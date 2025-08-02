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
            if (HEIGHT != WIDTH) @compileError("multiplicative identity matrix must be square");
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
    };
}
