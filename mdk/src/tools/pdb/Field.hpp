#pragma once
#include <string>
#include <vector>

namespace mdk::pdb {
    class Field {
    public:
        virtual void read(const std::string &s) = 0;
        virtual void write(std::string &s) const = 0;
    };

    class Integer: public Field {
    private:
        int i, j, offset;
        int *v;

    public:
        Integer(int i, int j, int& v, int offset = 0):
            i{i}, j{j}, offset{offset}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    class Real: public Field {
    private:
        int i, j, n, m;
        double *v, scalar;

    public:
        Real(int i, int j, int n, int m, double& v, double scalar = 1.0):
            i{i}, j{j}, n{n}, m{m}, v{&v}, scalar{scalar} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    enum Mode {
        Exact = 0x1, Trim = 0x2, Left = 0x4, Right = 0x8
    };

    class String: public Field {
    private:
        int i, j, n;
        std::string *v;
        int mode;

    public:
        String(int i, int j, std::string& v, int n = -1, int mode = Trim|Right):
            i{i}, j{j}, n{n}, v{&v}, mode{mode} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    class Char: public Field {
    private:
        int i;
        char *v;

    public:
        Char(int i, char& v):
            i{i}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    class SymOp: public Field {
    private:
        int i, j;
        std::string *v;

    public:
        SymOp(int i, int j, std::string& v):
            i{i}, j{j}, v{&v} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };

    class Literal: public Field {
    private:
        int i, j;
        std::string lit;

    public:
        Literal(int i, int j, std::string lit):
            i{i}, j{j}, lit{move(lit)} {};

        void read(const std::string &s) override;
        void write(std::string &s) const override;
    };
}