#include"Node.h"

Node::Node() {
	this->x = 0;
	this->y = 0;
	this->t = 0;
	this->BC = false;
}

void Node::set_x(double p) {
	this->x = p;
}

void Node::set_y(double p) {
	this->y = p;
}

void Node::set_BC(bool p) {
	this->BC = p;
}

void Node::set_t(double p) {
	this->t = p;
}

double Node::get_x() {
	return this->x;
}

double Node::get_y() {
	return this->y;
}

double Node::get_t() {
	return this->t;
}

bool Node::get_BC() {
	return this->BC;
}