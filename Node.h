#pragma once
class Node {
protected:
	double x, y;
	double t;
	bool BC;
public:
	Node();
	void set_x(double);
	double get_x();
	void set_y(double);
	double get_y();
	void set_BC(bool);
	bool get_BC();
	void set_t(double);
	double get_t();

};