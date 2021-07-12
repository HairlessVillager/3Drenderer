#include <cstdio>
#include <graphics.h>
#include <cmath>
#include <list>
#include <vector>
#include "3Drenderer.h"
using namespace std;

//#define NEGE
//#define DEBUG
#define my_iszero(x)	(fabs(x) < 1e-5)

const int SCREEN_WIDTH = 640.0;
const int SCREEN_HEIGHT = 480.0;
const int FPS = 60;

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
		Matrix4 rotateX(double theta) {
			Matrix4 res(1.0);
			res.value[1][1] = cos(theta);
			res.value[1][2] = sin(theta);
			res.value[2][1] = -sin(theta);
			res.value[2][2] = cos(theta);
			return *this * res;
		}
		Matrix4 rotateY(double theta) {
			Matrix4 res(1.0);
			res.value[0][0] = cos(theta);
			res.value[0][2] = -sin(theta);
			res.value[2][0] = sin(theta);
			res.value[2][2] = cos(theta);
			return *this * res;
		}
		Matrix4 rotateZ(double theta) {
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
			return Arrow3(begin.matrix4mul(mat), end.matrix4mul(mat));
		}
		~Arrow3() {
			;
		}
};

class Object {
	private :
		bool cut(Arrow3& arr, Camera& camera);
		bool cutNear(Arrow3& arr, double n);
		bool cutFar(Arrow3& arr, double f);
		bool cutLeft(Arrow3& arr, double fov, double aspect);
		bool cutRight(Arrow3& arr, double fov, double aspect);
		bool cutUp(Arrow3& arr, double fov);
		bool cutDown(Arrow3& arr, double fov);
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
			Arrow3 arr(Vector3(67,
			                   -103.67744862860489,
			                   156.31998191486267),
			           Vector3(66.9999999999999986,
			                   -19.153796280464981,
			                   -24.941575492467308));
			printf("\n%lf\t%lf\t%lf\n",
			       arr.begin.x,
			       arr.begin.y,
			       arr.begin.z);
			printf("%lf\t%lf\t%lf\n\n",
			       arr.end.x,
			       arr.end.y,
			       arr.end.z);
			cutUp(arr, 120.0 * M_PI / 180.0);
			printf("\n%lf\t%lf\t%lf\n",
			       arr.begin.x,
			       arr.begin.y,
			       arr.begin.z);
			printf("%lf\t%lf\t%lf\n\n",
			       arr.end.x,
			       arr.end.y,
			       arr.end.z);
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
			return product.matrix4mul(Matrix4(1.0).translate(
			                              Vector3(
			                                  - side * a * 0.5,
			                                  - side * b * 0.5,
			                                  0.0)));
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
		const double _sight = 300.0;
		const double _fov = 120.0 * M_PI / 180.0;
		const double _aspect = (double)SCREEN_WIDTH / SCREEN_HEIGHT;
		const double _near = 2.0;
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

bool Object::cutNear(Arrow3& arr, double n) {
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

bool Object::cutFar(Arrow3& arr, double f) {
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	if(begin.x <= f && end.x <= f) {
		return true;
	} else if(begin.x > f && end.x > f) {
		return false;
	} else {
		double x = f;
		double y = ((f - end.x) * begin.y - (f - begin.x) * end.y)
		           / (begin.x - end.x);
		double z = ((f - end.x) * begin.z - (f - begin.x) * end.z)
		           / (begin.x - end.x);
		if(begin.x < f) {
			arr.end.x = x;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = x;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

bool Object::cutUp(Arrow3& arr, double fov) {
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	double k_0 = atan(fov * 0.5);
	double k_1 = begin.z / begin.x;
	double k_2 = end.z / end.x;

	if(k_1 <= k_0 && k_2 <= k_0) {
		return true;
	} else if(k_1 > k_0 && k_2 > k_0) {
		return false;
	} else {
		double x, y, z;
		x = (begin.x * end.z - end.x * begin.z)
		    / (end.z - begin.z + k_0 * (begin.x - end.x));
		z = k_0 * x;
		if(my_iszero(begin.x - end.x)) {
			y = (end.y * (z - begin.z) - begin.y * (z - end.z))
			    / (end.z - begin.z);
		} else {
			y = (begin.y * (x - end.x) - end.y * (x - begin.x))
			    / (begin.x - end.x);
		}

		if(k_2 > k_0) {
			arr.end.x = x;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = x;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

bool Object::cutDown(Arrow3& arr, double fov) {
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	double k_0 = atan(- fov * 0.5);
	double k_1 = begin.z / begin.x;
	double k_2 = end.z / end.x;

	if(k_1 >= k_0 && k_2 >= k_0) {
		return true;
	} else if(k_1 < k_0 && k_2 < k_0) {
		return false;
	} else {
		double x, y, z;
		x = (begin.x * end.z - end.x * begin.z)
		    / (end.z - begin.z + k_0 * (begin.x - end.x));
		z = k_0 * x;
		if(my_iszero(begin.x - end.x)) {
			y = (end.y * (z - begin.z) - begin.y * (z - end.z))
			    / (end.z - begin.z);
		} else {
			y = (begin.y * (x - end.x) - end.y * (x - begin.x))
			    / (begin.x - end.x);
		}

		if(k_2 < k_0) {
			arr.end.x = x;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = x;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

bool Object::cutLeft(Arrow3& arr, double fov, double aspect) {
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	double k_0 = atan(fov * 0.5) * aspect;
	double k_1 = begin.y / begin.x;
	double k_2 = end.y / end.x;

	if(k_1 <= k_0 && k_2 <= k_0) {
		return true;
	} else if(k_1 > k_0 && k_2 > k_0) {
		return false;
	} else {
		double x, y, z;
		x = (begin.x * end.y - end.x * begin.y)
		    / (end.y - begin.y + k_0 * (begin.x - end.x));
		y = k_0 * x;
		if(my_iszero(begin.x - end.x)) {
			z = (begin.z * (y - end.y) - end.z * (y - begin.y))
			    / (begin.y - end.y);
		} else {
			z = (begin.z * (x - end.x) - end.z * (x - begin.x))
			    / (begin.x - end.x);
		}

		if(k_2 > k_0) {
			arr.end.x = x;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = x;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

bool Object::cutRight(Arrow3& arr, double fov, double aspect) {
	Vector3 begin = arr.begin;
	Vector3 end = arr.end;
	double k_0 = atan(- fov * 0.5) * aspect;
	double k_1 = begin.y / begin.x;
	double k_2 = end.y / end.x;

	if(k_1 >= k_0 && k_2 >= k_0) {
		return true;
	} else if(k_1 < k_0 && k_2 < k_0) {
		return false;
	} else {
		double x, y, z;
		x = (begin.x * end.y - end.x * begin.y)
		    / (end.y - begin.y + k_0 * (begin.x - end.x));
		y = k_0 * x;
		if(my_iszero(begin.x - end.x)) {
			z = (begin.z * (y - end.y) - end.z * (y - begin.y))
			    / (begin.y - end.y);
		} else {
			z = (begin.z * (x - end.x) - end.z * (x - begin.x))
			    / (begin.x - end.x);
		}

		if(k_2 < k_0) {
			arr.end.x = x;
			arr.end.y = y;
			arr.end.z = z;
		} else {
			arr.begin.x = x;
			arr.begin.y = y;
			arr.begin.z = z;
		}
		return true;
	}
}

bool Object::cut(Arrow3& arr, Camera& camera) {
	bool succ;
#ifdef DEBUG
	printf("before:\n%.2lf\t%.2lf\t%.2lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%.2lf\t%.2lf\t%.2lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#endif
	succ = cutNear(arr, camera.getNear());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutNear():\n%.2lf\t%.2lf\t%.2lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%.2lf\t%.2lf\t%.2lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#endif
	succ = cutFar(arr, camera.getFar());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutFar():\n%lf\t%lf\t%lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%lf\t%lf\t%lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#ifndef NEGE
	getch();
#endif
#endif
	succ = cutUp(arr, camera.getFov());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutUp():\n%lf\t%lf\t%lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%lf\t%lf\t%lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#ifndef NEGE
	getch();
#endif
#endif
	succ = cutDown(arr, camera.getFov());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutDown():\n%.2lf\t%.2lf\t%.2lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%.2lf\t%.2lf\t%.2lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#endif
	succ = cutLeft(arr, camera.getFov(), camera.getAspect());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutLeft():\n%.2lf\t%.2lf\t%.2lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%.2lf\t%.2lf\t%.2lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
#endif
	succ = cutRight(arr, camera.getFov(), camera.getAspect());
	if(!succ) {
		return false;
	}
#ifdef DEBUG
	printf("after cutRight():\n%.2lf\t%.2lf\t%.2lf\n",
	       arr.begin.x,
	       arr.begin.y,
	       arr.begin.z);
	printf("%.2lf\t%.2lf\t%.2lf\n\n",
	       arr.end.x,
	       arr.end.y,
	       arr.end.z);
	printf("==================\n");
#endif
	return true;
}

void Object::render(Camera& camera) {
	Object self_copy;
	Vector3 pos = camera.pos;
	Vector3 ang = camera.ang;

	// 世界坐标 -> 相机坐标
	self_copy = this->matrix4mul(Matrix4(1.0)
	                             .translate(-pos)
	                             .rotateZ(-ang.z)
	                             .rotateY(+ang.y)
	                             .rotateX(-ang.x));
	for(int i = 0; i < (int)self_copy.edges.size(); i++) {
		// 裁剪

		bool succ = cut(self_copy.edges[i], camera);

		if(succ) {
			// 相机坐标 -> 屏幕坐标，渲染线段
			Vector3 begin = self_copy.edges[i].begin;
			Vector3 end = self_copy.edges[i].end;
#ifndef NEGE
			line(
			    + begin.y * camera.getSight() / begin.x + SCREEN_WIDTH / 2.0,
			    - begin.z * camera.getSight() / begin.x + SCREEN_HEIGHT / 2.0,
			    + end.y * camera.getSight() / end.x + SCREEN_WIDTH / 2.0,
			    - end.z * camera.getSight() / end.x + SCREEN_HEIGHT / 2.0
			);
#ifdef DEBUG
			getch();
#endif
#endif
		}
	}
}

class Game {
	private :
		ObjectCreator _creator;
		Camera _camera;
		vector<Object> _obj_list;

		void interact(Camera& camera) {
			const double SPEED_0 = 1.0;
			const double SIGHT_SPEED_0 = M_PI / 180.0;
			double speed = SPEED_0;
			double sight_speed = SIGHT_SPEED_0;

			if(keystate(17)) {
				speed *= 2.0;
				sight_speed *= 2.0;
			}
			if(keystate(87)) {
				camera.pos.x += speed * cos(camera.ang.z);
				camera.pos.y += speed * sin(camera.ang.z);
			}
			if(keystate(83)) {
				camera.pos.x -= speed * cos(camera.ang.z);
				camera.pos.y -= speed * sin(camera.ang.z);
			}
			if(keystate(65)) {
				camera.pos.x += speed * cos(camera.ang.z - M_PI / 2);
				camera.pos.y += speed * sin(camera.ang.z - M_PI / 2);
			}
			if(keystate(68)) {
				camera.pos.x -= speed * cos(camera.ang.z - M_PI / 2);
				camera.pos.y -= speed * sin(camera.ang.z - M_PI / 2);
			}
			if(keystate(key_shift)) {
				camera.pos.z -= speed;
			} else if(keystate(key_space)) {
				camera.pos.z += speed;
			}

			if(camera.ang.y <= M_PI / 2 && keystate(key_up)) {
				camera.ang.y += sight_speed;
			} else if(camera.ang.y >= - M_PI / 2 && keystate(key_down)) {
				camera.ang.y -= sight_speed;
			}
			if(keystate(key_right)) {
				camera.ang.z += sight_speed;
				if(camera.ang.z > M_PI) {
					camera.ang.z -= 2 * M_PI;
				}
			} else if(keystate(key_left)) {
				camera.ang.z -= sight_speed ;
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
			Object ball, wall;

			ball = _creator.ball(10.0);
			wall = _creator.grid(10.0, 10, 10)
			       .matrix4mul(Matrix4(1.0).rotateY(M_PI * 0.5));

			for(; is_run(); delay_fps(FPS)) {
#ifndef EGE
				cleardevice();
				interact(_camera);
#endif
				for(int i = 0; i < (int)_obj_list.size(); i++) {
					_obj_list[i].render(_camera);
				}
				ball = ball.matrix4mul(Matrix4(1.0)
				                       .rotateZ(0.005));
				wall = wall.matrix4mul(Matrix4(1.0)
				                       .rotateX(0.005));
				ball.matrix4mul(Matrix4(1.0)
				                .translate(Vector3(40.0, 0.0, 30.0)))
				.render(_camera);
				wall.matrix4mul(Matrix4(1.0)
				                .rotateX(M_PI * 0.25)
				                .translate(Vector3(80.0, 0.0, 20.0)))
				.render(_camera);
#ifndef NEGE
				printDebugInfo();
#endif
			}
		}
	public :
		Game() {
#ifndef NEGE
			initgraph(SCREEN_WIDTH, SCREEN_HEIGHT,
			          INIT_TOPMOST | INIT_RENDERMANUAL);
#endif
			_obj_list.push_back(_creator.grid(10.0, 20, 20));
			_camera.pos = Vector3(0.0, 0.0, 5.0);
		}
		void run() {
			mainLoop();
		}
		void debug() {
			Object obj;

			_camera.pos = Vector3(85.49, 42.10, 67.00);
			_camera.ang = Vector3(0.0 / 180.0 * M_PI,
			                      -90.0 / 180.0 * M_PI,
			                      -155.0 / 180.0 * M_PI);

			obj = _creator.grid(10.0, 10, 10)
			      .matrix4mul(Matrix4(1.0)
			                  .rotateY(M_PI * 0.5)
			                  .rotateX(M_PI * 0.25)
			                  .translate(Vector3(80.0, 0.0, 20.0)));
#ifndef NEGE
			cleardevice();
#endif
			for(int i = 0; i < (int)_obj_list.size(); i++) {
				_obj_list[i].render(_camera);
			}
			//obj.render(_camera);
#ifndef NEGE
			printDebugInfo();
			getch();
#endif
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
