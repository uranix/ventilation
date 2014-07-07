#ifndef __VEC_H__
#define __VEC_H__

#include <cmath>

namespace dir {
	enum Direction {
		X = 0,
		Y = 1,
		Z = 2,
		DIR_BEGIN = X,
		DIR_END = 3,
	};

	enum DirectionRange {
		DIRECTIONS
	};

	class DirectionIterator {
		Direction dir;
	public:
		DirectionIterator(Direction dir) : dir(dir) { }
		DirectionIterator &operator++() {
			dir = static_cast<Direction>(static_cast<unsigned>(dir) + 1);
			return *this;
		}
		Direction operator*() const {
			return dir;
		}
		bool operator!=(DirectionIterator other) {
			return this->dir != other.dir;
		}
	};

	inline DirectionIterator begin(DirectionRange) { return DirectionIterator(DIR_BEGIN); }
	inline DirectionIterator end(DirectionRange) { return DirectionIterator(DIR_END); }

	inline char to_char(Direction dir) {
		if (dir == dir::X)
			return 'X';
		if (dir == dir::Y)
			return 'Y';
		return 'Z';
	}

	inline Direction from_char(char cdir) {
		if (cdir == 'x' || cdir == 'X')
			return dir::X;
		if (cdir == 'y' || cdir == 'Y')
			return dir::Y;
		return dir::Z;
	}

	template<typename T>
	T &select(Direction dir, T &i, T &j, T &k) {
		if (dir == dir::X)
			return i;
		if (dir == dir::Y)
			return j;
		return k;
	}

	template<typename T>
	const T &select(Direction dir, const T &i, const T &j, const T &k) {
		if (dir == dir::X)
			return i;
		if (dir == dir::Y)
			return j;
		return k;
	}
};

struct vec {
    double x; //!< vector x component
    double y; //!< vector y component
    double z; //!< vector z component
    /** Construct zero vector */
    vec() : x(0), y(0), z(0) { }
    /** Construct a vector with all components equal to a */
    explicit vec(double a) : x(a), y(a), z(a) { }
    /** Construct vector from another vector r */
    vec(const vec &r) : x(r.x), y(r.y), z(r.z) { }
    /** Construct vector from components */
    vec(double x, double y, double z) : x(x), y(y), z(z) { }

    double &operator()(dir::Direction d) {
        if (d == dir::X)
            return x;
        if (d == dir::Y)
            return y;
        return z;
    }
    const double &operator()(dir::Direction d) const {
        if (d == dir::X)
            return x;
        if (d == dir::Y)
            return y;
        return z;
    }

#define BINOP(op) \
    const vec operator op(const vec &r) const { \
        return vec(x op r.x, y op r.y, z op r.z); \
    }

    /** Add two vectors */
    BINOP(+)
    /** Subtract two vectors */
    BINOP(-)
    /** Element-wise multiply two vectors */
    BINOP(*)
    /** Element-wise divide two vectors */
    BINOP(/)

#undef BINOP

    /** Multiply vector by number a */
    const vec operator *(const double a) const {
        return vec(a * x, a * y, a * z);
    }
    /** Divide vector by number a */
    const vec operator /(const double a) const {
        return vec(x / a, y / a, z / a);
    }
    /** Negate a vector */
    const vec operator -() const {
        return vec(-x, -y, -z);
    }
    /** Scale vector by number a */
    vec &operator *=(double a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }
    /** Cross-product vector with another vector b */
    const vec operator %(const vec &b) const {
        const vec &a = *this;
        return vec(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
    }
    /** Return square of vector's euclidean norm */
    double norm2() const {
        return x * x + y * y + z * z;
    }
    /** Return vector's euclidean norm */
    double norm() const {
        return sqrt(norm2());
    }
    /** Dot-product vector with another vector b */
    double dot(const vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }

#define COMPOP(op) \
    vec &operator op(const vec &b) { \
        x op b.x; \
        y op b.y; \
        z op b.z; \
        return *this; \
    }

    /** Add another vector to this vector */
    COMPOP(+=)
        /** Subtract another vector from this vector */
    COMPOP(-=)
        /** Element-wise multiply this vector by another vector */
    COMPOP(*=)
        /** Element-wise divide this vector by another vector */
    COMPOP(/=)

#undef COMPOP
};

/** Multiply vector b by number a */
inline const vec operator *(double a, const vec &b) {
    return b * a;
}

/** Return dot-product of two vectors */
inline double dot(const vec &a, const vec &b) {
    return a.dot(b);
}

/** Return triple product of three vectors, i.e. a.dot(b.cross(c)) */
inline double triple(const vec &a, const vec &b, const vec &c) {
    return a.dot(b % c);
}

/** Return squared norm of vector a */
inline double norm2(const vec &a) {
    return a.norm2();
}

/** Return norm of vector a */
inline double norm(const vec &a) {
    return a.norm();
}

#include <ostream>
/** Pretty print vector r to output stream o */
inline std::ostream &operator <<(std::ostream &o, const vec &r) {
    return o << "(" << r.x << ", " << r.y << ", " << r.z << ")";
}

namespace dir {
	inline vec to_vec(Direction dir) {
		if (dir == dir::X)
			return vec(1, 0, 0);
		if (dir == dir::Y)
			return vec(0, 1, 0);
		return vec(0, 0, 1);
	}
}

#endif
