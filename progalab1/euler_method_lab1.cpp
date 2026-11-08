#include <cmath>
#include <deque>
#include <vector>
#include <iostream>
#include "fsolve.hpp"

const double M1 = 1.0;
const double M2 = 1.0;
const double L1 = 1.0;
const double L2 = 1.0;
const double G = 9.80665;
const double pi = 3.141592653589793;

double dt = 0.0005;

struct State
{
	double x1, vx1, x2, vx2, y1, vy1, y2, vy2;
};

struct T_vect {
	double T1, T2;
};

T_vect calc_forces_from_coords(const State& state) {
	double a1 = -state.x1 * state.x1 / (M1 * L1) - state.y1 * state.y1 / (M1 * L1);
	double b1 = (state.x2 - state.x1) * state.x1 / (M1 * L2) + (state.y2 - state.y1) * state.y1 / (M1 * L2);

	double a2 = (state.x1 * (state.x2 - state.x1) + state.y1 * (state.y2 - state.y1)) / (M1 * L1);
	double b2 = (1.0 / M1 + 1.0 / M2) / L2 * ((state.y2 - state.y1) * (state.y2 - state.y1) + (state.x2 - state.x1) * (state.x2 - state.x1));

	double c1 = G * state.y1 - state.vy1 * state.vy1 - state.vx1 * state.vx1;
	double c2 = -(state.vx2 - state.vx1) * (state.vx2 - state.vx1) - (state.vy2 - state.vy1) * (state.vy2 - state.vy1);

	double T2 = (c2 - c1 * a2 / a1) / (b2 - b1 * a2 / a1);
	double T1 = (c1 - b1 * T2) / a1;

	return { T1, T2 };
}

State calc_coords_from_forces(const T_vect& T, const State& prev_state) {
	double a_1 = (T.T1 / L1 + T.T2 / L2) * dt / M1;
	double b_1 = -dt * T.T2 / L2 / M1;
	double c_1 = 1.;
	double d_1 = 0.;

	double a_2 = 1.;
	double b_2 = 0.;
	double c_2 = -dt;
	double d_2 = 0.;

	double a_3 = -dt * T.T2 / L2 / M2;
	double b_3 = dt * T.T2 / L2 / M2;
	double c_3 = 0.;
	double d_3 = 1.;

	double a_4 = 0.;
	double b_4 = 1.;
	double c_4 = 0.;
	double d_4 = -dt;

	double h_1 = prev_state.vx1;
	double h_2 = prev_state.x1;
	double h_3 = prev_state.vx2;
	double h_4 = prev_state.x2;

	double x1 = (-b_4 * c_3 * d_2 * h_1 + b_3 * c_4 * d_2 * h_1 + b_4 * c_2 * d_3 * h_1 - b_2 * c_4 * d_3 * h_1 - b_3 * c_2 * d_4 * h_1 + b_2 * c_3 * d_4 * h_1 + b_4 * c_3 * d_1 * h_2 - b_3 * c_4 * d_1 * h_2 - b_4 * c_1 * d_3 * h_2 + b_1 * c_4 * d_3 * h_2 + b_3 * c_1 * d_4 * h_2 - b_1 * c_3 * d_4 * h_2 - b_4 * c_2 * d_1 * h_3 + b_2 * c_4 * d_1 * h_3 + b_4 * c_1 * d_2 * h_3 - b_1 * c_4 * d_2 * h_3 - b_2 * c_1 * d_4 * h_3 + b_1 * c_2 * d_4 * h_3 + b_3 * c_2 * d_1 * h_4 - b_2 * c_3 * d_1 * h_4 - b_3 * c_1 * d_2 * h_4 + b_1 * c_3 * d_2 * h_4 + b_2 * c_1 * d_3 * h_4 - b_1 * c_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);
	double x2 = (a_4 * c_3 * d_2 * h_1 - a_3 * c_4 * d_2 * h_1 - a_4 * c_2 * d_3 * h_1 + a_2 * c_4 * d_3 * h_1 + a_3 * c_2 * d_4 * h_1 - a_2 * c_3 * d_4 * h_1 - a_4 * c_3 * d_1 * h_2 + a_3 * c_4 * d_1 * h_2 + a_4 * c_1 * d_3 * h_2 - a_1 * c_4 * d_3 * h_2 - a_3 * c_1 * d_4 * h_2 + a_1 * c_3 * d_4 * h_2 + a_4 * c_2 * d_1 * h_3 - a_2 * c_4 * d_1 * h_3 - a_4 * c_1 * d_2 * h_3 + a_1 * c_4 * d_2 * h_3 + a_2 * c_1 * d_4 * h_3 - a_1 * c_2 * d_4 * h_3 - a_3 * c_2 * d_1 * h_4 + a_2 * c_3 * d_1 * h_4 + a_3 * c_1 * d_2 * h_4 - a_1 * c_3 * d_2 * h_4 - a_2 * c_1 * d_3 * h_4 + a_1 * c_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);

	double vx1 = (-a_4 * b_3 * d_2 * h_1 + a_3 * b_4 * d_2 * h_1 + a_4 * b_2 * d_3 * h_1 - a_2 * b_4 * d_3 * h_1 - a_3 * b_2 * d_4 * h_1 + a_2 * b_3 * d_4 * h_1 + a_4 * b_3 * d_1 * h_2 - a_3 * b_4 * d_1 * h_2 - a_4 * b_1 * d_3 * h_2 + a_1 * b_4 * d_3 * h_2 + a_3 * b_1 * d_4 * h_2 - a_1 * b_3 * d_4 * h_2 - a_4 * b_2 * d_1 * h_3 + a_2 * b_4 * d_1 * h_3 + a_4 * b_1 * d_2 * h_3 - a_1 * b_4 * d_2 * h_3 - a_2 * b_1 * d_4 * h_3 + a_1 * b_2 * d_4 * h_3 + a_3 * b_2 * d_1 * h_4 - a_2 * b_3 * d_1 * h_4 - a_3 * b_1 * d_2 * h_4 + a_1 * b_3 * d_2 * h_4 + a_2 * b_1 * d_3 * h_4 - a_1 * b_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);
	double vx2 = (a_4 * b_3 * c_2 * h_1 - a_3 * b_4 * c_2 * h_1 - a_4 * b_2 * c_3 * h_1 + a_2 * b_4 * c_3 * h_1 + a_3 * b_2 * c_4 * h_1 - a_2 * b_3 * c_4 * h_1 - a_4 * b_3 * c_1 * h_2 + a_3 * b_4 * c_1 * h_2 + a_4 * b_1 * c_3 * h_2 - a_1 * b_4 * c_3 * h_2 - a_3 * b_1 * c_4 * h_2 + a_1 * b_3 * c_4 * h_2 + a_4 * b_2 * c_1 * h_3 - a_2 * b_4 * c_1 * h_3 - a_4 * b_1 * c_2 * h_3 + a_1 * b_4 * c_2 * h_3 + a_2 * b_1 * c_4 * h_3 - a_1 * b_2 * c_4 * h_3 - a_3 * b_2 * c_1 * h_4 + a_2 * b_3 * c_1 * h_4 + a_3 * b_1 * c_2 * h_4 - a_1 * b_3 * c_2 * h_4 - a_2 * b_1 * c_3 * h_4 + a_1 * b_2 * c_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);


	h_1 = prev_state.vy1 - dt * G;
	h_2 = prev_state.y1;
	h_3 = prev_state.vy2 - dt * G;
	h_4 = prev_state.y2;

	double y1 = (-b_4 * c_3 * d_2 * h_1 + b_3 * c_4 * d_2 * h_1 + b_4 * c_2 * d_3 * h_1 - b_2 * c_4 * d_3 * h_1 - b_3 * c_2 * d_4 * h_1 + b_2 * c_3 * d_4 * h_1 + b_4 * c_3 * d_1 * h_2 - b_3 * c_4 * d_1 * h_2 - b_4 * c_1 * d_3 * h_2 + b_1 * c_4 * d_3 * h_2 + b_3 * c_1 * d_4 * h_2 - b_1 * c_3 * d_4 * h_2 - b_4 * c_2 * d_1 * h_3 + b_2 * c_4 * d_1 * h_3 + b_4 * c_1 * d_2 * h_3 - b_1 * c_4 * d_2 * h_3 - b_2 * c_1 * d_4 * h_3 + b_1 * c_2 * d_4 * h_3 + b_3 * c_2 * d_1 * h_4 - b_2 * c_3 * d_1 * h_4 - b_3 * c_1 * d_2 * h_4 + b_1 * c_3 * d_2 * h_4 + b_2 * c_1 * d_3 * h_4 - b_1 * c_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);
	double y2 = (a_4 * c_3 * d_2 * h_1 - a_3 * c_4 * d_2 * h_1 - a_4 * c_2 * d_3 * h_1 + a_2 * c_4 * d_3 * h_1 + a_3 * c_2 * d_4 * h_1 - a_2 * c_3 * d_4 * h_1 - a_4 * c_3 * d_1 * h_2 + a_3 * c_4 * d_1 * h_2 + a_4 * c_1 * d_3 * h_2 - a_1 * c_4 * d_3 * h_2 - a_3 * c_1 * d_4 * h_2 + a_1 * c_3 * d_4 * h_2 + a_4 * c_2 * d_1 * h_3 - a_2 * c_4 * d_1 * h_3 - a_4 * c_1 * d_2 * h_3 + a_1 * c_4 * d_2 * h_3 + a_2 * c_1 * d_4 * h_3 - a_1 * c_2 * d_4 * h_3 - a_3 * c_2 * d_1 * h_4 + a_2 * c_3 * d_1 * h_4 + a_3 * c_1 * d_2 * h_4 - a_1 * c_3 * d_2 * h_4 - a_2 * c_1 * d_3 * h_4 + a_1 * c_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);
	double vy1 = (-a_4 * b_3 * d_2 * h_1 + a_3 * b_4 * d_2 * h_1 + a_4 * b_2 * d_3 * h_1 - a_2 * b_4 * d_3 * h_1 - a_3 * b_2 * d_4 * h_1 + a_2 * b_3 * d_4 * h_1 + a_4 * b_3 * d_1 * h_2 - a_3 * b_4 * d_1 * h_2 - a_4 * b_1 * d_3 * h_2 + a_1 * b_4 * d_3 * h_2 + a_3 * b_1 * d_4 * h_2 - a_1 * b_3 * d_4 * h_2 - a_4 * b_2 * d_1 * h_3 + a_2 * b_4 * d_1 * h_3 + a_4 * b_1 * d_2 * h_3 - a_1 * b_4 * d_2 * h_3 - a_2 * b_1 * d_4 * h_3 + a_1 * b_2 * d_4 * h_3 + a_3 * b_2 * d_1 * h_4 - a_2 * b_3 * d_1 * h_4 - a_3 * b_1 * d_2 * h_4 + a_1 * b_3 * d_2 * h_4 + a_2 * b_1 * d_3 * h_4 - a_1 * b_2 * d_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);


	double vy2 = (a_4 * b_3 * c_2 * h_1 - a_3 * b_4 * c_2 * h_1 - a_4 * b_2 * c_3 * h_1 + a_2 * b_4 * c_3 * h_1 + a_3 * b_2 * c_4 * h_1 - a_2 * b_3 * c_4 * h_1 - a_4 * b_3 * c_1 * h_2 + a_3 * b_4 * c_1 * h_2 + a_4 * b_1 * c_3 * h_2 - a_1 * b_4 * c_3 * h_2 - a_3 * b_1 * c_4 * h_2 + a_1 * b_3 * c_4 * h_2 + a_4 * b_2 * c_1 * h_3 - a_2 * b_4 * c_1 * h_3 - a_4 * b_1 * c_2 * h_3 + a_1 * b_4 * c_2 * h_3 + a_2 * b_1 * c_4 * h_3 - a_1 * b_2 * c_4 * h_3 - a_3 * b_2 * c_1 * h_4 + a_2 * b_3 * c_1 * h_4 + a_3 * b_1 * c_2 * h_4 - a_1 * b_3 * c_2 * h_4 - a_2 * b_1 * c_3 * h_4 + a_1 * b_2 * c_3 * h_4) / (a_4 * b_3 * c_2 * d_1 - a_3 * b_4 * c_2 * d_1 - a_4 * b_2 * c_3 * d_1 + a_2 * b_4 * c_3 * d_1 + a_3 * b_2 * c_4 * d_1 - a_2 * b_3 * c_4 * d_1 - a_4 * b_3 * c_1 * d_2 + a_3 * b_4 * c_1 * d_2 + a_4 * b_1 * c_3 * d_2 - a_1 * b_4 * c_3 * d_2 - a_3 * b_1 * c_4 * d_2 + a_1 * b_3 * c_4 * d_2 + a_4 * b_2 * c_1 * d_3 - a_2 * b_4 * c_1 * d_3 - a_4 * b_1 * c_2 * d_3 + a_1 * b_4 * c_2 * d_3 + a_2 * b_1 * c_4 * d_3 - a_1 * b_2 * c_4 * d_3 - a_3 * b_2 * c_1 * d_4 + a_2 * b_3 * c_1 * d_4 + a_3 * b_1 * c_2 * d_4 - a_1 * b_3 * c_2 * d_4 - a_2 * b_1 * c_3 * d_4 + a_1 * b_2 * c_3 * d_4);

	return { x1, vx1, x2, vx2, y1, vy1, y2, vy2 };
}

double t_stop = 10.0;
double angle1 = pi / 2;
double angle2 = pi / 2;
int steps = static_cast<int>(t_stop / dt);

std::vector<double> t(steps);

State state0 = { L1 * sin(angle1), 0.0, L1 * sin(angle1) + L2 * sin(angle2), 0.0,
				-L1 * cos(angle1), 0.0, -L1 * cos(angle1) - L2 * cos(angle2), 0.0 };
T_vect T_0 = calc_forces_from_coords(state0);

std::vector<State> y = { state0 };
int i;

void fcn(int n, double* x, double* fvec) {
	T_vect T = { x[0], x[1] };
	State inner_state1 = calc_coords_from_forces(T, y[i - 1]);

	fvec[0] = inner_state1.x1 * inner_state1.x1 + inner_state1.y1 * inner_state1.y1 - L1 * L1;
	fvec[1] = (inner_state1.x2 - inner_state1.x1) * (inner_state1.x2 - inner_state1.x1) + (inner_state1.y2 - inner_state1.y1) * (inner_state1.y2 - inner_state1.y1) - L2 * L2;

}

int main() {
	i = 0;
	std::cout << dt * i << " " << y[i].x1 << " " << y[i].y1 << " " << y[i].x2 << " " << y[i].y2 << " " << y[i].vx1 << " " << y[i].vy1 << " " << y[i].vx2 << " " << y[i].vy2 << "\n";
	for (i = 1; i < steps; i++) {
		State current_state = calc_coords_from_forces(T_0, y[i - 1]);
		T_vect T_1 = calc_forces_from_coords(current_state);

		double diff = sqrt(pow(T_0.T1 - T_1.T1, 2) + pow(T_0.T2 - T_1.T2, 2));
		if (diff > 1e-14) {
			double x[2] = { T_0.T1, T_0.T2 };
			double tol = 1e-14;
			int n = 2;
			double fvec[2];
			double wa[19];

			int status = fsolve(fcn, n, x, fvec, tol, wa, 19);

			T_0 = T_vect{ x[0], x[1] };
		}

		y.push_back(calc_coords_from_forces(T_0, y[i - 1]));
		std::cout << dt * i << " " << y[i].x1 << " " << y[i].y1 << " " << y[i].x2 << " " << y[i].y2 << " " << y[i].vx1 << " " << y[i].vy1 << " " << y[i].vx2 << " " << y[i].vy2 << "\n";
		T_0 = calc_forces_from_coords(y[i]);
	}

	return 0;
}