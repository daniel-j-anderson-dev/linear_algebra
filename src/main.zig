const std = @import("std");

const linear_algebra = @import("linear_algebra");
const Matrix = linear_algebra.Matrix;

pub fn main() !void {
    var general_purpose_allocator = std.heap.GeneralPurposeAllocator(.{}).init;
    defer _ = general_purpose_allocator.deinit();
    const allocator = general_purpose_allocator.allocator();
    _ = allocator;
}
