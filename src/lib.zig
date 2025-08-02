const builtin = @import("builtin");

pub const matrix = @import("matrix.zig");
pub const Matrix = matrix.Matrix;

pub const IN_DEBUG_MODE = builtin.mode == .Debug;
