#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "MPC.h"
#include "json.hpp"
#include "utils.h"

#define MPH2MPS 0.44704
// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}



int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"]; v *= MPH2MPS ; 
          double current_steering_angle = j[1]["steering_angle"];
          double current_throttle = j[1]["throttle"];
          

          /********************************************************************
          * Calculate steering angle and throttle using MPC.                  *
          * Both are in between [-1, 1].                                      *
          ********************************************************************/
          
          // Convert reference trajectory points to car's coordinate system.
          // This will make life 10x-100x easier but it won't guarantee 
          // that another Ariane5(https://around.com/ariane.html) will not happen.
          Eigen::VectorXd ref_xs(ptsx.size());
          Eigen::VectorXd ref_ys(ptsy.size());
          for(size_t i=0; i<ptsx.size(); i++){
            ref_xs(i) = (ptsx[i] - px) * cos(-psi) - (ptsy[i] - py) * sin(-psi) ;
            ref_ys(i) = (ptsx[i] - px) * sin(-psi) + (ptsy[i] - py) * cos(-psi) ;
          }  
          
          // Fit a polynomial to the reference trajectory.
          Eigen::VectorXd reference_curve = polyfit(ref_xs,ref_ys,3);

          // Find current cross track error CTE. 
          double ct_err = polyeval(reference_curve,0);
          
          // Find current heading error which is the difference of tangent 
          // direction at current position and the current heading, which is
          // 0 in the car's frame of reference.
          double ref_heading = atan(find_slope(reference_curve,0));           
          double psi_err = 0 - ref_heading;
          
          cout << "cte= " << ct_err << " psi_err= " << psi_err << " ";

          // Adjust for latency of 100 ms.
          const double latency = 0.1;
          const double Lf = 2.67;
          const double x_adjusted = 0 + v*latency;
          const double y_adjusted = 0;
          const double psi_adjusted = 0 - v * current_steering_angle * latency / Lf;
          const double v_adjusted = v + current_throttle * latency;
          const double ct_err_adjusted = ct_err + v * sin(psi_err) * latency;
          const double psi_err_adjusted = psi_err + psi_adjusted; 

          // Fill up the state vector.
          Eigen::VectorXd state(6);
          state << x_adjusted, y_adjusted, psi_adjusted, v_adjusted, ct_err_adjusted, psi_err_adjusted;

          // Call the solver and get the steering angle, throttle that minimize the 
          // the deviation from the reference state over the entire prediction horizon.
          auto vars = mpc.Solve(state,reference_curve);
          
          double steer_value = -vars[0]/deg2rad(25);
          double throttle_value = vars[1];

          cout << "steer_value= " << steer_value << " throttle= " << throttle_value << endl <<endl;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          for(size_t i=2; i<vars.size(); i+=2){
             mpc_x_vals.push_back(vars[i]);
             mpc_y_vals.push_back(vars[i+1]);
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for(int i=0; i<ref_xs.size(); i++){
            next_x_vals.push_back(ref_xs(i));
            next_y_vals.push_back(ref_ys(i));
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
