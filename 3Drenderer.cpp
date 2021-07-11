#include <cstdio>
#include <graphics.h>
#include <cmath>
#include <list>
#include <vector>
#include "3Drenderer.h"
using namespace std;

const int SCREEN_WIDTH = 640.0;
const int SCREEN_HEIGHT = 480.0;

class Vector3 {
	public :
		double x;
		double y;
		double z;
		Vector3() {
			x = y = z = 0.0;
		}
		Vector3(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
		}
		Vector3 operator+(const Vector3& v) {
			return Vector3(this->x + v.x, this->y + v.y, this->z + v.z);
		}
		Vector3 operator-(const Vector3& v) {
			return Vector3(this->x - v.x, this->y - v.y, this->z - v.z);
		}
		Vector3 operator-() {
			return Vector3(-x, -y, -z);
		}
		double operator*(const Vector3& v) {
			return this->x * v.x + this->y * v.y + this->z * v.z;
		}
		Vector3 matrix4mul(const Matrix4& mat);
		double length() {
			return sqrt((this->x) * (this->x)
			            + (this->y) * (this->y)
			            + (this->z) * (this->z));
		}
		double length2() {
			return (this->x) * (this->x)
			       + (this->y) * (this->y)
			       + (this->z) * (this->z);
		}
		~Vector3() {
			;
		}
};

class Matrix4 {
	public :
		double value[4][4];
		Matrix4(double x) {
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					if(i == j) {
						value[i][j] = x;
					} else {
						value[i][j] = 0.0;
					}
				}
			}
		}
		Matrix4 operator*(const Matrix4& mat) {
			Matrix4 res(0.0);
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					for(int k = 0; k < 4; k++) {
						res.value[i][j] += this->value[i][k] * mat.value[k][j];
					}
				}
			}
			return res;
		}
		Matrix4 translate(Vector3 dv) {
			Matrix4 res(1.0);
			res.value[3][0] = dv.x;
			res.value[3][1] = dv.y;
			res.value[3][2] = dv.z;
			return *this * res;
		}
		Matrix4 scale(Vector3 sv) {
			Matrix4 res(1.0);
			res.value[0][0] = sv.x;
			res.value[1][1] = sv.y;
			res.value[2][2] = sv.z;
			return *this * res;
		}
		Matrix4 revolveX(double theta) {
			Matrix4 res(1.0);
			res.value[1][1] = cos(theta);
			res.value[1][2] = sin(theta);
			res.value[2][1] = -sin(theta);
			res.value[2][2] = cos(theta);
			return *this * res;
		}
		Matrix4 revolveY(double theta) {
			Matrix4 res(1.0);
			res.value[0][0] = cos(theta);
			res.value[0][2] = -sin(theta);
			res.value[2][0] = sin(theta);
			res.value[2][2] = cos(theta);
			return *this * res;
		}
		Matrix4 revolveZ(double theta) {
			Matrix4 res(1.0);
			res.value[0][0] = cos(theta);
			res.value[0][1] = sin(theta);
			res.value[1][0] = -sin(theta);
			res.value[1][1] = cos(theta);
			return *this * res;
		}
		void debug() {
			for(int j = 0; j < 4; j++) {
				for(int i = 0; i < 4; i++) {
					printf("%.1lf\t", value[i][j]);
				}
				printf("\n");
			}
		}
		~Matrix4() {
			;
		}
};

Vector3 Vector3::matrix4mul(const Matrix4& mat) {
	double vector4[1][4] = {{x, y, z, 1.0}};
	double res[1][4] = {{0.0, 0.0, 0.0, 0.0}};
	for(int i = 0; i < 1; i++) {
		for(int j = 0; j < 4; j++) {
			for(int k = 0; k < 4; k++) {
				res[i][j] += vector4[i][k] * mat.value[k][j];
			}
		}
	}
	return Vector3(res[0][0], res[0][1], res[0][2]);
}

class Arrow3 {
	public :
		Vector3 begin;
		Vector3 end;
		Arrow3() {
			begin.x = begin.y = begin.z = 0.0;
			end.x = end.y = end.z = 0.0;
		}
		Arrow3(Vector3 begin, Vector3 end) {
			this->begin = begin;
			this->end = end;
		}
		Arrow3 matrix4mul(const Matrix4& mat) {
			begin = begin.matrix4mul(mat);
			end = end.matrix4mul(mat);
			return Arrow3(begin, end);
		}
		~Arrow3() {
			;
		}
};

class Object {
	private :
		bool cut(Arrow3& arr, double n);
		vector<Arrow3> edges;

	public :
		Object() {
			;
		}
		void render(Camera& camera);
		void addEdge(Arrow3 edge) {
			edges.push_back((edge));
		}
		Object operator+(const Object& obj) {
			Object product;
			product.edges.resize(this->edges.size() + obj.edges.size());
			int i = 0;
			int j = 0;
			int k = 0;
			while(j < (int)this->edges.size()) {
				product.edges[i++] = this->edges[j++];
			}
			while(k < (int)obj.edges.size()) {
				product.edges[i++] = obj.edges[k++];
			}
			return product;
		}
		Object matrix4mul(const Matrix4& mat) {
			Object product;
			for(int i = 0; i < (int)this->edges.size(); i++) {
				product.addEdge(this->edges[i].matrix4mul(mat));
			}
			return product;
		}
		void debug() {
			printf("x\ty\tz\n");
			for(int i = 0; i < (int)edges.size(); i++) {
				printf("index = %d :\n", i);
				printf("%.1lf\t%.1lf\t%.1lf\n",
				       edges[i].begin.x,
				       edges[i].begin.y,
				       edges[i].begin.z);
				printf("%.1lf\t%.1lf\t%.1lf\n",
				       edges[i].end.x,
				       edges[i].end.y,
				       edges[i].end.z);
			}
		}
		~Object() {
			;
		}
};

class ObjectCreator {
	public :
		ObjectCreator() {
			;
		}
		Object cube(double side) {
			Object product;
			double half = side * 0.5;
			// 顶部
			for(int i = 0; i < 4; i++) {
				Vector3 begin((i & 02) ? -half : +half,
				              (i & 02) ? -half : +half,
				              +half);
				Vector3 end((i & 01) ? +half : -half,
				            (i & 01) ? -half : +half,
				            +half);
				product.addEdge(Arrow3(begin, end));
			}
			// 腰部
			for(int i = 0; i < 4; i++) {
				Vector3 begin((i & 01) ? +half : -half,
				              (i & 02) ? +half : -half,
				              +half);
				Vector3 end((i & 01) ? +half : -half,
				            (i & 02) ? +half : -half,
				            -half);
				product.addEdge(Arrow3(begin, end));
			}
			// 底部
			for(int i = 0; i < 4; i++) {
				Vector3 begin((i & 02) ? -half : +half,
				              (i & 02) ? -half : +half,
				              -half);
				Vector3 end((i & 01) ? +half : -half,
				            (i & 01) ? -half : +half,
				            -half);
				product.addEdge(Arrow3(begin, end));
			}
			return product;
		}
		Object grid(double side, int a, int b) {
			Object product;
			for(int yi = 0; yi <= b; yi++) {
				Vector3 begin(0.0, yi * side, 0.0);
				Vector3 end(a * side, yi * side, 0.0);
				product.addEdge(Arrow3(begin, end));
			}
			for(int xi = 0; xi <= a; xi++) {
				Vector3 begin(xi * side, 0.0, 0.0);
				Vector3 end(xi * side, b * side, 0.0);
				product.addEdge(Arrow3(begin, end));
			}
			product.matrix4mul(Matrix4(1.0).translate(
			                       Vector3(
			                           - side * a * 0.5,
			                           - side * b * 0.5,
			                           0.0)));
			return product;
		}
		Object ball(double radius) {
			Object product;
			double pi_10 = M_PI * 0.1;
			double pi_8 = M_PI * 0.125;
			// 添加经线 20
			for(int i = 0; i < 20; i++) {
				for(int j = 0; j < 4; j++) {
					double x_0 = radius * cos(i * pi_10);
					double y_0 = radius * sin(i * pi_10);
					Vector3 begin(x_0 * cos(j * pi_8),
					              y_0 * cos(j * pi_8),
					              radius * sin(j * pi_8));
					Vector3 end(x_0 * cos((j + 1) * pi_8),
					            y_0 * cos((j + 1) * pi_8),
					            radius * sin((j + 1) * pi_8));
					product.addEdge(Arrow3(begin, end));
				}
				for(int j = 0; j > -4; j--) {
					double x_0 = radius * cos(i * pi_10);
					double y_0 = radius * sin(i * pi_10);
					Vector3 begin(x_0 * cos(j * pi_8),
					              y_0 * cos(j * pi_8),
					              radius * sin(j * pi_8));
					Vector3 end(x_0 * cos((j - 1) * pi_8),
					            y_0 * cos((j - 1) * pi_8),
					            radius * sin((j - 1) * pi_8));
					product.addEdge(Arrow3(begin, end));
				}
			}
			// 添加纬线 7
			for(int j = -3; j <= 3; j++) {
				double z_0 = radius * sin(j * pi_8);
				for(int i = 0; i < 20; i++) {
					Vector3 begin(radius * cos(i * pi_10) * cos(j * pi_8),
					              radius * sin(i * pi_10) * cos(j * pi_8),
					              z_0);
					Vector3 end(radius * cos((i + 1) * pi_10) * cos(j * pi_8),
					            radius * sin((i + 1) * pi_10) * cos(j * pi_8),
					            z_0);
					product.addEdge(Arrow3(begin, end));
				}
			}
			return product;
		}
		~ObjectCreator() {
			;
		}
};

class Camera {
	private :
		const double _sight = 150.0;
		const double _fov = 50.0;
		const double _aspect = 1.0;
		const double _near = 0.1;
		const double _far = 200.0;

	public :
		Vector3 pos;
		Vector3 ang;

		Camera() {
			pos.x = pos.y = pos.z = 0.0;
			ang.x = ang.y = ang.z = 0.0;
		}

		const double getSight() {
			return this->_sight;
		}
		const double getFov() {
			return this->_fov;
		}
		const double getAspect() {
			return this->_aspect;
		}
		const double getNear() {
			return this->_near;
		}
		const double getFar() {
			return this->_far;
		}
		~Camera() {
			;
		}
};

bool Object::cut(Arrow3& arr, double n) {
	// todo : 这里先裁 x <= near 的部分，其余的以后再说
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	if(begin.x > n && end.x > n) {
		return true;
	} else if(begin.x <= n && end.x <= n) {
		return false;
	} else {
		double y = ((n - end.x) * begin.y - (n - begin.x) * end.y)
		           / (begin.x - end.x);
		double z = ((n - end.x) * begin.z - (n - begin.x) * end.z)
		           / (begin.x - end.x);
		if(begin.x > n) {
			arr.end.x = n;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = n;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

void Object::render(Camera& camera) {
	Object self_copy = *this;
	Vector3 pos = camera.pos;
	Vector3 ang = camera.ang;
	double sight = camera.getSight();

	// 世界坐标 -> 相机坐标
	self_copy.matrix4mul(Matrix4(1.0)
	                     .translate(-pos)
	                     .revolveZ(-ang.z)
	                     .revolveY(+ang.y)
	                     .revolveX(-ang.x));
	for(int i = 0; i < (int)self_copy.edges.size(); i++) {
		// 裁剪
		bool succ = cut(self_copy.edges[i], camera.getNear());
		if(succ) {
			// 相机坐标 -> 屏幕坐标，渲染线段
			Vector3 begin = self_copy.edges[i].begin;
			Vector3 end = self_copy.edges[i].end;
			line(
			    + begin.y * sight / begin.x + SCREEN_WIDTH / 2.0,
			    - begin.z * sight / begin.x + SCREEN_HEIGHT / 2.0,
			    + end.y * sight / end.x + SCREEN_WIDTH / 2.0,
			    - end.z * sight / end.x + SCREEN_HEIGHT / 2.0
			);
		}
	}
}

class Game {
	private :
		ObjectCreator _creator;
		Camera _camera;
		vector<Object> _obj_list;
		const int FPS = 60;

		void interact(Camera& camera) {
			const double SPEED = 1.0;
			const double SIGHT_SPEED = M_PI / 180.0;
			if(keystate(87)) {
				camera.pos.x += SPEED * cos(camera.ang.z);
				camera.pos.y += SPEED * sin(camera.ang.z);
			}
			if(keystate(83)) {
				camera.pos.x -= SPEED * cos(camera.ang.z);
				camera.pos.y -= SPEED * sin(camera.ang.z);
			}
			if(keystate(65)) {
				camera.pos.x += SPEED * cos(camera.ang.z - M_PI / 2);
				camera.pos.y += SPEED * sin(camera.ang.z - M_PI / 2);
			}
			if(keystate(68)) {
				camera.pos.x -= SPEED * cos(camera.ang.z - M_PI / 2);
				camera.pos.y -= SPEED * sin(camera.ang.z - M_PI / 2);
			}
			if(keystate(key_shift)) {
				camera.pos.z -= SPEED ;
			} else if(keystate(key_space)) {
				camera.pos.z += SPEED ;
			}

			if(camera.ang.y <= M_PI / 2 && keystate(key_up)) {
				camera.ang.y += SIGHT_SPEED ;
			} else if(camera.ang.y >= - M_PI / 2 && keystate(key_down)) {
				camera.ang.y -= SIGHT_SPEED ;
			}
			if(keystate(key_right)) {
				camera.ang.z += SIGHT_SPEED ;
				if(camera.ang.z > M_PI) {
					camera.ang.z -= 2 * M_PI;
				}
			} else if(keystate(key_left)) {
				camera.ang.z -= SIGHT_SPEED ;
				if(camera.ang.z < - M_PI) {
					camera.ang.z += 2 * M_PI;
				}
			}
		}
		void printDebugInfo() {
			xyprintf(0, 0, "%.2f", getfps());
			xyprintf(0, 20, "%.2lf %.2lf %.2lf",
			         _camera.pos.x,
			         _camera.pos.y,
			         _camera.pos.z);
			xyprintf(0, 40, "%.2lf %.2lf %.2lf",
			         _camera.ang.x / M_PI * 180,
			         _camera.ang.y / M_PI * 180,
			         _camera.ang.z / M_PI * 180);
		}
		void mainLoop() {
			Object obj;
			obj = _creator.ball((20.0)) + _creator.ball((10.0));
			obj = obj.matrix4mul(Matrix4(1.0).translate(
			                         Vector3(0.0, 0.0, 30.0)));
			for(; is_run(); delay_fps(FPS)) {
				cleardevice();
				interact(_camera);
				for(int i = 0; i < (int)_obj_list.size(); i++) {
					_obj_list[i].render(_camera);
				}
				obj = obj.matrix4mul(Matrix4(1.0).revolveZ(0.005));
				obj.render(_camera);
				printDebugInfo();
			}
		}
	public :
		Game() {
			initgraph(SCREEN_WIDTH, SCREEN_HEIGHT, INIT_RENDERMANUAL);
			_obj_list.push_back(_creator.grid(10.0, 20, 20));
			_camera.pos = Vector3(0.0, 0.0, 5.0);
		}
		void run() {
			mainLoop();
		}
		~Game() {
			;
		}
};

int main() {
	Game game;
	game.run();
	return 0;
}
