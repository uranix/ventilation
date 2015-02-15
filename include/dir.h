#ifndef __DIR_H__
#define __DIR_H__

#include <cmath>

namespace dir {
    enum Direction {
        X = 0,
        Y = 1,
        Z = 2,
        DIR_BEGIN = X,
        DIR_END = 3,
    };

    enum Side {
        BEG = 0,
        END = 1,
        SIDE_BEGIN = BEG,
        SIDE_END = 2
    };

    inline Side flip(Side v) {
        if (v == BEG)
            return END;
        return BEG;
    }

    enum DirectionRange {
        DIRECTIONS
    };

    enum SideRange {
        SIDES
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

    class SideIterator {
        Side s;
    public:
        SideIterator(Side s) : s(s) { }
        SideIterator &operator++() {
            s = static_cast<Side>(static_cast<unsigned>(s) + 1);
            return *this;
        }
        Side operator*() const {
            return s;
        }
        bool operator!=(SideIterator other) {
            return this->s != other.s;
        }
    };

    inline DirectionIterator begin(DirectionRange) { return DirectionIterator(DIR_BEGIN); }
    inline DirectionIterator end(DirectionRange) { return DirectionIterator(DIR_END); }

    inline SideIterator begin(SideRange) { return SideIterator(SIDE_BEGIN); }
    inline SideIterator end(SideRange) { return SideIterator(SIDE_END); }

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
}

#endif
