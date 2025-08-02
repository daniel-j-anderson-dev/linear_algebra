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
    };
}
