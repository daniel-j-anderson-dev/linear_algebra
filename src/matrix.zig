const builtin = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;

const linear_algebra = @import("linear_algebra");
const IN_DEBUG_MODE = linear_algebra.IN_DEBUG_MODE;

pub fn Matrix(T: type, H: usize, W: usize) type {
    if (H == 0 or W == 0) @compileError("0 size Matrix is not supported");

    return struct {
        // fields
        rows: [HEIGHT][WIDTH]T,

        // constants
        pub const HEIGHT = H;
        pub const WIDTH = W;
        pub const SIZE = HEIGHT * WIDTH;
        pub const IS_SQUARE = HEIGHT == WIDTH;

        // types
        const Self = @This();

        // helpers
        pub fn outOfBounds(row: usize, column: usize) bool {
            return row >= HEIGHT or column >= WIDTH;
        }
        pub fn assertInBounds(row: usize, column: usize) void {
            if (outOfBounds(row, column)) @panic("index out of bounds");
        }

        // constructors
        pub const ZEROS = Self{ .rows = [1][WIDTH]T{[_]T{0} ** WIDTH} ** HEIGHT };
        pub fn multiplicativeIdentity() Self {
            if (!IS_SQUARE) @compileError("multiplicative identity matrix must be square");

            var id = Self.ZEROS;
            for (0..HEIGHT) |i| {
                const id_element = id.getMutableUnchecked(i, i);
                id_element.* = 1;
            }
            return id;
        }
        pub fn fromArray(rows: [HEIGHT][WIDTH]T) Self {
            return Self{ .rows = rows };
        }
        pub fn clone(self: *const Self) Self {
            return Self.fromArray(&self.rows);
        }
        pub fn matrixOfCofactors(self: *const Self) Self {
            var c = Self.ZEROS;
            for (0..HEIGHT) |i| {
                for (0..WIDTH) |j| {
                    const cofactor_element = c.getMutableUnchecked(i, j);
                    cofactor_element.* = self.cofactorUnchecked(i, j);
                }
            }
            return c;
        }

        // getters
        pub fn height(_: *const Self) usize {
            return HEIGHT;
        }
        pub fn width(_: *const Self) usize {
            return WIDTH;
        }
        pub fn get(self: *const Self, row: usize, column: usize) ?*const T {
            if (outOfBounds(row, column)) return null;
            return &self.rows[row][column];
        }
        pub fn getMutable(self: *const Self, row: usize, column: usize) ?*T {
            if (outOfBounds(row, column)) return null;
            return &self.rows[row][column];
        }
        pub fn getUnchecked(self: *const Self, row: usize, column: usize) *const T {
            if (IN_DEBUG_MODE) assertInBounds(row, column);
            return &self.rows[row][column];
        }
        pub fn getMutableUnchecked(self: *const Self, row: usize, column: usize) *T {
            if (IN_DEBUG_MODE) assertInBounds(row, column);
            return &self.rows[row][column];
        }

        // operators
        pub fn transpose(self: *const Matrix(T, HEIGHT, WIDTH)) Matrix(T, WIDTH, HEIGHT) {
            var transposed = Matrix(T, WIDTH, HEIGHT).ZEROS;
            for (0..HEIGHT) |i| {
                for (0..WIDTH) |j| {
                    const self_element = self.getUnchecked(i, j);
                    const transpose_element = transposed.getMutableUnchecked(j, i);

                    transpose_element.* = self_element.*;
                }
            }
            return transposed;
        }
        pub fn add(lhs: *const Self, rhs: *const Self) Self {
            var sum = Self.ZEROS;
            for (0..HEIGHT) |i| {
                for (0..WIDTH) |j| {
                    const lhs_element = lhs.getUnchecked(i, j);
                    const rhs_element = rhs.getUnchecked(i, j);
                    const sum_element = sum.getMutableUnchecked(i, j);

                    sum_element.* = lhs_element.* + rhs_element.*;
                }
            }
            return sum;
        }
        pub fn negate(self: *const Self) Self {
            var negated = self.clone();
            for (0..HEIGHT) |i| {
                for (0..WIDTH) |j| {
                    const self_element = self.getUnchecked(i, j);
                    const negated_element = negated.getMutUnchecked(i, j);
                    negated_element.* = -(self_element.*);
                }
            }
            return negated;
        }
        pub fn subtract(lhs: *const Self, rhs: *const Self) Self {
            return lhs.add(rhs.negate());
        }
        pub fn multiply(
            comptime RHS_WIDTH: usize,
            lhs: *const Matrix(T, HEIGHT, WIDTH),
            rhs: *const Matrix(T, WIDTH, RHS_WIDTH),
        ) Matrix(T, HEIGHT, RHS_WIDTH) {
            var product = Matrix(T, HEIGHT, RHS_WIDTH).ZEROS;
            for (0..HEIGHT) |lhs_row_index| {
                for (0..RHS_WIDTH) |rhs_column_index| {
                    const dot_product = product.getMutableUnchecked(lhs_row_index, rhs_column_index);
                    for (0..WIDTH) |i| {
                        const lhs_element = lhs.getUnchecked(lhs_row_index, i);
                        const rhs_element = rhs.getUnchecked(i, rhs_column_index);

                        dot_product.* += lhs_element.* * rhs_element.*;
                    }
                }
            }
            return product;
        }
        pub fn scalar_multiply(self: *const Self, scalar: *const T) Self {
            var product = Self.ZEROS;
            for (0..HEIGHT) |i| {
                for (0..WIDTH) |j| {
                    const product_element = product.getMutableUnchecked(i, j);
                    const self_element = self.getUnchecked(i, j);

                    product_element.* = scalar * self_element.*;
                }
            }
            return product;
        }
        pub fn determinant(self: *const Self) T {
            if (!IS_SQUARE) @compileError("Determinant is only defined for square matrices");

            if (WIDTH == 2 and HEIGHT == 2) {
                const a = self.getUnchecked(0, 0).*;
                const b = self.getUnchecked(0, 1).*;
                const c = self.getUnchecked(1, 0).*;
                const d = self.getUnchecked(1, 1).*;
                return a * d - c * b;
            }
            if (WIDTH == 1 and HEIGHT == 1) return self.getUnchecked(0, 0).*;

            var sum: T = 0;
            const i = 0;
            for (0..WIDTH) |j| {
                const self_element = self.getUnchecked(i, j);
                const self_cofactor = self.cofactor(i, j);

                sum = sum + (self_element.* * self_cofactor);
            }
            return sum;
        }
        pub fn minor(self: *const Self, row: usize, column: usize) ?Matrix(T, HEIGHT - 1, WIDTH - 1) {
            if (outOfBounds(row, column)) return null;
            return self.minorUnchecked(row, column);
        }
        pub fn minorUnchecked(self: *const Self, row: usize, column: usize) Matrix(T, HEIGHT - 1, WIDTH - 1) {
            if (IN_DEBUG_MODE) assertInBounds(row, column);

            var minor_ = Matrix(T, HEIGHT - 1, WIDTH - 1).ZEROS;

            var minor_row: usize = 0;
            for (0..HEIGHT) |i| {
                if (i == row) continue;

                var minor_column: usize = 0;
                for (0..WIDTH) |j| {
                    if (j == column) continue;

                    const minor_element = minor_.getMutableUnchecked(minor_row, minor_column);
                    const self_element = self.getUnchecked(i, j);

                    minor_element.* = self_element.*;

                    minor_column += 1;
                }
                minor_row += 1;
            }
            return minor_;
        }
        pub fn cofactor(self: *const Self, row: usize, column: usize) ?T {
            if (outOfBounds(row, column)) return null;
            return self.cofactorUnchecked(row, column);
        }
        pub fn cofactorUnchecked(self: *const Self, row: usize, column: usize) T {
            const sign: T = if ((row + column) % 2 == 0) 1 else -1;
            return sign * self.minorUnchecked(row, column).determinant();
        }
        pub fn inverse(self: *const Self) ?Self {
            const det = self.determinant();
            if (det == 0) return null;
            const C = self.matrixOfCofactors();
            return C.transpose().scalar_multiply(det);
        }
        pub fn inverseUnchecked(self: *const Self) Self {
            const det = self.determinant();
            if (IN_DEBUG_MODE) if (det == 0) @panic("Inverse does not exists for matrices with determinant of zero");
            const C = self.matrixOfCofactors();
            return C.transpose().scalar_multiply(det);
        }
    };
}
