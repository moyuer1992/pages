"use strict";

var _get = function get(object, property, receiver) { if (object === null) object = Function.prototype; var desc = Object.getOwnPropertyDescriptor(object, property); if (desc === undefined) { var parent = Object.getPrototypeOf(object); if (parent === null) { return undefined; } else { return get(parent, property, receiver); } } else if ("value" in desc) { return desc.value; } else { var getter = desc.get; if (getter === undefined) { return undefined; } return getter.call(receiver); } };

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

//////////////////////////////////////////////////////////////////////////////
//
//  Angel.js
//
//////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
//
//  Helper functions
//

function _argumentsToArray(args) {
  return [].concat.apply([], Array.prototype.slice.apply(args));
}

//----------------------------------------------------------------------------

function radians(degrees) {
  return degrees * Math.PI / 180.0;
}

//----------------------------------------------------------------------------
//
//  Vector Constructors
//

function vec2() {
  var result = _argumentsToArray(arguments);

  switch (result.length) {
    case 0:
      result.push(0.0);
    case 1:
      result.push(0.0);
  }

  return result.splice(0, 2);
}

function vec3() {
  var result = _argumentsToArray(arguments);

  switch (result.length) {
    case 0:
      result.push(0.0);
    case 1:
      result.push(0.0);
    case 2:
      result.push(0.0);
  }

  return result.splice(0, 3);
}

function vec4() {
  var result = _argumentsToArray(arguments);

  switch (result.length) {
    case 0:
      result.push(0.0);
    case 1:
      result.push(0.0);
    case 2:
      result.push(0.0);
    case 3:
      result.push(1.0);
  }

  return result.splice(0, 4);
}

//----------------------------------------------------------------------------
//
//  Matrix Constructors
//

function mat2() {
  var v = _argumentsToArray(arguments);

  var m = [];
  switch (v.length) {
    case 0:
      v[0] = 1;
    case 1:
      m = [vec2(v[0], 0.0), vec2(0.0, v[0])];
      break;

    default:
      m.push(vec2(v));v.splice(0, 2);
      m.push(vec2(v));
      break;
  }

  m.matrix = true;

  return m;
}

//----------------------------------------------------------------------------

function mat3() {
  var v = _argumentsToArray(arguments);

  var m = [];
  switch (v.length) {
    case 0:
      v[0] = 1;
    case 1:
      m = [vec3(v[0], 0.0, 0.0), vec3(0.0, v[0], 0.0), vec3(0.0, 0.0, v[0])];
      break;

    default:
      m.push(vec3(v));v.splice(0, 3);
      m.push(vec3(v));v.splice(0, 3);
      m.push(vec3(v));
      break;
  }

  m.matrix = true;

  return m;
}

//----------------------------------------------------------------------------

function mat4() {
  var v = _argumentsToArray(arguments);

  var m = [];
  switch (v.length) {
    case 0:
      v[0] = 1;
    case 1:
      m = [vec4(v[0], 0.0, 0.0, 0.0), vec4(0.0, v[0], 0.0, 0.0), vec4(0.0, 0.0, v[0], 0.0), vec4(0.0, 0.0, 0.0, v[0])];
      break;

    default:
      m.push(vec4(v));v.splice(0, 4);
      m.push(vec4(v));v.splice(0, 4);
      m.push(vec4(v));v.splice(0, 4);
      m.push(vec4(v));
      break;
  }

  m.matrix = true;

  return m;
}

//----------------------------------------------------------------------------
//
//  Generic Mathematical Operations for Vectors and Matrices
//

function equal(u, v) {
  if (u.length != v.length) {
    return false;
  }

  if (u.matrix && v.matrix) {
    for (var i = 0; i < u.length; ++i) {
      if (u[i].length != v[i].length) {
        return false;
      }
      for (var j = 0; j < u[i].length; ++j) {
        if (u[i][j] !== v[i][j]) {
          return false;
        }
      }
    }
  } else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
    return false;
  } else {
    for (var i = 0; i < u.length; ++i) {
      if (u[i] !== v[i]) {
        return false;
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------

function add(u, v) {
  var result = [];

  if (u.matrix && v.matrix) {
    if (u.length != v.length) {
      throw "add(): trying to add matrices of different dimensions";
    }

    for (var i = 0; i < u.length; ++i) {
      if (u[i].length != v[i].length) {
        throw "add(): trying to add matrices of different dimensions";
      }
      result.push([]);
      for (var j = 0; j < u[i].length; ++j) {
        result[i].push(u[i][j] + v[i][j]);
      }
    }

    result.matrix = true;

    return result;
  } else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
    throw "add(): trying to add matrix and non-matrix variables";
  } else {
    if (u.length != v.length) {
      throw "add(): vectors are not the same dimension";
    }

    for (var i = 0; i < u.length; ++i) {
      result.push(u[i] + v[i]);
    }

    return result;
  }
}

//----------------------------------------------------------------------------

function subtract(u, v) {
  var result = [];

  if (u.matrix && v.matrix) {
    if (u.length != v.length) {
      throw "subtract(): trying to subtract matrices" + " of different dimensions";
    }

    for (var i = 0; i < u.length; ++i) {
      if (u[i].length != v[i].length) {
        throw "subtract(): trying to subtact matrices" + " of different dimensions";
      }
      result.push([]);
      for (var j = 0; j < u[i].length; ++j) {
        result[i].push(u[i][j] - v[i][j]);
      }
    }

    result.matrix = true;

    return result;
  } else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
    throw "subtact(): trying to subtact  matrix and non-matrix variables";
  } else {
    if (u.length != v.length) {
      throw "subtract(): vectors are not the same length";
    }

    for (var i = 0; i < u.length; ++i) {
      result.push(u[i] - v[i]);
    }

    return result;
  }
}

//----------------------------------------------------------------------------

function mult(u, v) {
  var result = [];

  if (u.matrix && v.matrix) {
    if (u.length != v.length) {
      throw "mult(): trying to add matrices of different dimensions";
    }

    for (var i = 0; i < u.length; ++i) {
      if (u[i].length != v[i].length) {
        throw "mult(): trying to add matrices of different dimensions";
      }
    }

    for (var i = 0; i < u.length; ++i) {
      result.push([]);

      for (var j = 0; j < v.length; ++j) {
        var sum = 0.0;
        for (var k = 0; k < u.length; ++k) {
          sum += u[i][k] * v[k][j];
        }
        result[i].push(sum);
      }
    }

    result.matrix = true;

    return result;
  }

  if (u.matrix && u.length == v.length) {
    for (var i = 0; i < v.length; i++) {
      var sum = 0.0;
      for (var j = 0; j < v.length; j++) {
        sum += u[i][j] * v[j];
      }
      result.push(sum);
    }
    return result;
  } else {
    if (u.length != v.length) {
      throw "mult(): vectors are not the same dimension";
    }

    for (var i = 0; i < u.length; ++i) {
      result.push(u[i] * v[i]);
    }

    return result;
  }
}

//----------------------------------------------------------------------------
//
//  Basic Transformation Matrix Generators
//

function translate(x, y, z) {
  if (Array.isArray(x) && x.length == 3) {
    z = x[2];
    y = x[1];
    x = x[0];
  }

  var result = mat4();
  result[0][3] = x;
  result[1][3] = y;
  result[2][3] = z;

  return result;
}

//----------------------------------------------------------------------------

function rotate(angle, axis) {
  if (!Array.isArray(axis)) {
    axis = [arguments[1], arguments[2], arguments[3]];
  }

  var v = normalize(axis);

  var x = v[0];
  var y = v[1];
  var z = v[2];

  var c = Math.cos(radians(angle));
  var omc = 1.0 - c;
  var s = Math.sin(radians(angle));

  var result = mat4(vec4(x * x * omc + c, x * y * omc - z * s, x * z * omc + y * s, 0.0), vec4(x * y * omc + z * s, y * y * omc + c, y * z * omc - x * s, 0.0), vec4(x * z * omc - y * s, y * z * omc + x * s, z * z * omc + c, 0.0), vec4());

  return result;
}

function rotateX(theta) {
  var c = Math.cos(radians(theta));
  var s = Math.sin(radians(theta));
  var rx = mat4(1.0, 0.0, 0.0, 0.0, 0.0, c, s, 0.0, 0.0, -s, c, 0.0, 0.0, 0.0, 0.0, 1.0);
  return rx;
}
function rotateY(theta) {
  var c = Math.cos(radians(theta));
  var s = Math.sin(radians(theta));
  var ry = mat4(c, 0.0, -s, 0.0, 0.0, 1.0, 0.0, 0.0, s, 0.0, c, 0.0, 0.0, 0.0, 0.0, 1.0);
  return ry;
}
function rotateZ(theta) {
  var c = Math.cos(radians(theta));
  var s = Math.sin(radians(theta));
  var rz = mat4(c, s, 0.0, 0.0, -s, c, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  return rz;
}

//----------------------------------------------------------------------------

function scalem(x, y, z) {
  if (Array.isArray(x) && x.length == 3) {
    z = x[2];
    y = x[1];
    x = x[0];
  }

  var result = mat4();
  result[0][0] = x;
  result[1][1] = y;
  result[2][2] = z;

  return result;
}

//----------------------------------------------------------------------------
//
//  ModelView Matrix Generators
//

function lookAt(eye, at, up) {
  if (!Array.isArray(eye) || eye.length != 3) {
    throw "lookAt(): first parameter [eye] must be an a vec3";
  }

  if (!Array.isArray(at) || at.length != 3) {
    throw "lookAt(): first parameter [at] must be an a vec3";
  }

  if (!Array.isArray(up) || up.length != 3) {
    throw "lookAt(): first parameter [up] must be an a vec3";
  }

  if (equal(eye, at)) {
    return mat4();
  }

  var v = normalize(subtract(at, eye)); // view direction vector
  var n = normalize(cross(v, up)); // perpendicular vector
  var u = normalize(cross(n, v)); // "new" up vector

  v = negate(v);

  var result = mat4(vec4(n, -dot(n, eye)), vec4(u, -dot(u, eye)), vec4(v, -dot(v, eye)), vec4());

  return result;
}

//----------------------------------------------------------------------------
//
//  Projection Matrix Generators
//

function ortho(left, right, bottom, top, near, far) {
  if (left == right) {
    throw "ortho(): left and right are equal";
  }
  if (bottom == top) {
    throw "ortho(): bottom and top are equal";
  }
  if (near == far) {
    throw "ortho(): near and far are equal";
  }

  var w = right - left;
  var h = top - bottom;
  var d = far - near;

  var result = mat4();
  result[0][0] = 2.0 / w;
  result[1][1] = 2.0 / h;
  result[2][2] = -2.0 / d;
  result[0][3] = -(left + right) / w;
  result[1][3] = -(top + bottom) / h;
  result[2][3] = -(near + far) / d;

  return result;
}

//----------------------------------------------------------------------------

function perspective(fovy, aspect, near, far) {
  var f = 1.0 / Math.tan(radians(fovy) / 2);
  var d = far - near;

  var result = mat4();
  result[0][0] = f / aspect;
  result[1][1] = f;
  result[2][2] = -(near + far) / d;
  result[2][3] = -2 * near * far / d;
  result[3][2] = -1;
  result[3][3] = 0.0;

  return result;
}

//----------------------------------------------------------------------------
//
//  Matrix Functions
//

function transpose(m) {
  if (!m.matrix) {
    return "transpose(): trying to transpose a non-matrix";
  }

  var result = [];
  for (var i = 0; i < m.length; ++i) {
    result.push([]);
    for (var j = 0; j < m[i].length; ++j) {
      result[i].push(m[j][i]);
    }
  }

  result.matrix = true;

  return result;
}

//----------------------------------------------------------------------------
//
//  Vector Functions
//

function dot(u, v) {
  if (u.length != v.length) {
    throw "dot(): vectors are not the same dimension";
  }

  var sum = 0.0;
  for (var i = 0; i < u.length; ++i) {
    sum += u[i] * v[i];
  }

  return sum;
}

//----------------------------------------------------------------------------

function negate(u) {
  var result = [];
  for (var i = 0; i < u.length; ++i) {
    result.push(-u[i]);
  }

  return result;
}

//----------------------------------------------------------------------------

function cross(u, v) {
  if (!Array.isArray(u) || u.length < 3) {
    throw "cross(): first argument is not a vector of at least 3";
  }

  if (!Array.isArray(v) || v.length < 3) {
    throw "cross(): second argument is not a vector of at least 3";
  }

  var result = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];

  return result;
}

//----------------------------------------------------------------------------

function length(u) {
  return Math.sqrt(dot(u, u));
}

//----------------------------------------------------------------------------

function normalize(u, excludeLastComponent) {
  if (excludeLastComponent) {
    var last = u.pop();
  }

  var len = length(u);

  if (!isFinite(len)) {
    throw "normalize: vector " + u + " has zero length";
  }

  for (var i = 0; i < u.length; ++i) {
    u[i] /= len;
  }

  if (excludeLastComponent) {
    u.push(last);
  }

  return u;
}

//----------------------------------------------------------------------------

function mix(u, v, s) {
  if (typeof s !== "number") {
    throw "mix: the last paramter " + s + " must be a number";
  }

  if (u.length != v.length) {
    throw "vector dimension mismatch";
  }

  var result = [];
  for (var i = 0; i < u.length; ++i) {
    result.push((1.0 - s) * u[i] + s * v[i]);
  }

  return result;
}

//----------------------------------------------------------------------------
//
// Vector and Matrix functions
//

function scale(s, u) {
  if (!Array.isArray(u)) {
    throw "scale: second parameter " + u + " is not a vector";
  }

  var result = [];
  for (var i = 0; i < u.length; ++i) {
    result.push(s * u[i]);
  }

  return result;
}

//----------------------------------------------------------------------------
//
//
//

function flatten(v) {
  if (v.matrix === true) {
    v = transpose(v);
  }

  var n = v.length;
  var elemsAreArrays = false;

  if (Array.isArray(v[0])) {
    elemsAreArrays = true;
    n *= v[0].length;
  }

  var floats = new Float32Array(n);

  if (elemsAreArrays) {
    var idx = 0;
    for (var i = 0; i < v.length; ++i) {
      for (var j = 0; j < v[i].length; ++j) {
        floats[idx++] = v[i][j];
      }
    }
  } else {
    for (var i = 0; i < v.length; ++i) {
      floats[i] = v[i];
    }
  }

  return floats;
}

//----------------------------------------------------------------------------

var sizeof = {
  'vec2': new Float32Array(flatten(vec2())).byteLength,
  'vec3': new Float32Array(flatten(vec3())).byteLength,
  'vec4': new Float32Array(flatten(vec4())).byteLength,
  'mat2': new Float32Array(flatten(mat2())).byteLength,
  'mat3': new Float32Array(flatten(mat3())).byteLength,
  'mat4': new Float32Array(flatten(mat4())).byteLength
};

// new functions 5/2/2015

// printing

function printm(m) {
  if (m.length == 2) for (var i = 0; i < m.length; i++) {
    console.log(m[i][0], m[i][1]);
  } else if (m.length == 3) for (var i = 0; i < m.length; i++) {
    console.log(m[i][0], m[i][1], m[i][2]);
  } else if (m.length == 4) for (var i = 0; i < m.length; i++) {
    console.log(m[i][0], m[i][1], m[i][2], m[i][3]);
  }
}
// determinants

function det2(m) {

  return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

function det3(m) {
  var d = m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[2][1] * m[1][0] - m[2][0] * m[1][1] * m[0][2] - m[1][0] * m[0][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1];
  return d;
}

function det4(m) {
  var m0 = [vec3(m[1][1], m[1][2], m[1][3]), vec3(m[2][1], m[2][2], m[2][3]), vec3(m[3][1], m[3][2], m[3][3])];
  var m1 = [vec3(m[1][0], m[1][2], m[1][3]), vec3(m[2][0], m[2][2], m[2][3]), vec3(m[3][0], m[3][2], m[3][3])];
  var m2 = [vec3(m[1][0], m[1][1], m[1][3]), vec3(m[2][0], m[2][1], m[2][3]), vec3(m[3][0], m[3][1], m[3][3])];
  var m3 = [vec3(m[1][0], m[1][1], m[1][2]), vec3(m[2][0], m[2][1], m[2][2]), vec3(m[3][0], m[3][1], m[3][2])];
  return m[0][0] * det3(m0) - m[0][1] * det3(m1) + m[0][2] * det3(m2) - m[0][3] * det3(m3);
}

function det(m) {
  if (m.matrix != true) console.log("not a matrix");
  if (m.length == 2) return det2(m);
  if (m.length == 3) return det3(m);
  if (m.length == 4) return det4(m);
}

//---------------------------------------------------------

// inverses

function inverse2(m) {
  var a = mat2();
  var d = det2(m);
  a[0][0] = m[1][1] / d;
  a[0][1] = -m[0][1] / d;
  a[1][0] = -m[1][0] / d;
  a[1][1] = m[0][0] / d;
  a.matrix = true;
  return a;
}

function inverse3(m) {
  var a = mat3();
  var d = det3(m);

  var a00 = [vec2(m[1][1], m[1][2]), vec2(m[2][1], m[2][2])];
  var a01 = [vec2(m[1][0], m[1][2]), vec2(m[2][0], m[2][2])];
  var a02 = [vec2(m[1][0], m[1][1]), vec2(m[2][0], m[2][1])];
  var a10 = [vec2(m[0][1], m[0][2]), vec2(m[2][1], m[2][2])];
  var a11 = [vec2(m[0][0], m[0][2]), vec2(m[2][0], m[2][2])];
  var a12 = [vec2(m[0][0], m[0][1]), vec2(m[2][0], m[2][1])];
  var a20 = [vec2(m[0][1], m[0][2]), vec2(m[1][1], m[1][2])];
  var a21 = [vec2(m[0][0], m[0][2]), vec2(m[1][0], m[1][2])];
  var a22 = [vec2(m[0][0], m[0][1]), vec2(m[1][0], m[1][1])];

  a[0][0] = det2(a00) / d;
  a[0][1] = -det2(a10) / d;
  a[0][2] = det2(a20) / d;
  a[1][0] = -det2(a01) / d;
  a[1][1] = det2(a11) / d;
  a[1][2] = -det2(a21) / d;
  a[2][0] = det2(a02) / d;
  a[2][1] = -det2(a12) / d;
  a[2][2] = det2(a22) / d;

  return a;
}

function inverse4(m) {
  var a = mat4();
  var d = det4(m);

  var a00 = [vec3(m[1][1], m[1][2], m[1][3]), vec3(m[2][1], m[2][2], m[2][3]), vec3(m[3][1], m[3][2], m[3][3])];
  var a01 = [vec3(m[1][0], m[1][2], m[1][3]), vec3(m[2][0], m[2][2], m[2][3]), vec3(m[3][0], m[3][2], m[3][3])];
  var a02 = [vec3(m[1][0], m[1][1], m[1][3]), vec3(m[2][0], m[2][1], m[2][3]), vec3(m[3][0], m[3][1], m[3][3])];
  var a03 = [vec3(m[1][0], m[1][1], m[1][2]), vec3(m[2][0], m[2][1], m[2][2]), vec3(m[3][0], m[3][1], m[3][2])];
  var a10 = [vec3(m[0][1], m[0][2], m[0][3]), vec3(m[2][1], m[2][2], m[2][3]), vec3(m[3][1], m[3][2], m[3][3])];
  var a11 = [vec3(m[0][0], m[0][2], m[0][3]), vec3(m[2][0], m[2][2], m[2][3]), vec3(m[3][0], m[3][2], m[3][3])];
  var a12 = [vec3(m[0][0], m[0][1], m[0][3]), vec3(m[2][0], m[2][1], m[2][3]), vec3(m[3][0], m[3][1], m[3][3])];
  var a13 = [vec3(m[0][0], m[0][1], m[0][2]), vec3(m[2][0], m[2][1], m[2][2]), vec3(m[3][0], m[3][1], m[3][2])];
  var a20 = [vec3(m[0][1], m[0][2], m[0][3]), vec3(m[1][1], m[1][2], m[1][3]), vec3(m[3][1], m[3][2], m[3][3])];
  var a21 = [vec3(m[0][0], m[0][2], m[0][3]), vec3(m[1][0], m[1][2], m[1][3]), vec3(m[3][0], m[3][2], m[3][3])];
  var a22 = [vec3(m[0][0], m[0][1], m[0][3]), vec3(m[1][0], m[1][1], m[1][3]), vec3(m[3][0], m[3][1], m[3][3])];
  var a23 = [vec3(m[0][0], m[0][1], m[0][2]), vec3(m[1][0], m[1][1], m[1][2]), vec3(m[3][0], m[3][1], m[3][2])];

  var a30 = [vec3(m[0][1], m[0][2], m[0][3]), vec3(m[1][1], m[1][2], m[1][3]), vec3(m[2][1], m[2][2], m[2][3])];
  var a31 = [vec3(m[0][0], m[0][2], m[0][3]), vec3(m[1][0], m[1][2], m[1][3]), vec3(m[2][0], m[2][2], m[2][3])];
  var a32 = [vec3(m[0][0], m[0][1], m[0][3]), vec3(m[1][0], m[1][1], m[1][3]), vec3(m[2][0], m[2][1], m[2][3])];
  var a33 = [vec3(m[0][0], m[0][1], m[0][2]), vec3(m[1][0], m[1][1], m[1][2]), vec3(m[2][0], m[2][1], m[2][2])];

  a[0][0] = det3(a00) / d;
  a[0][1] = -det3(a10) / d;
  a[0][2] = det3(a20) / d;
  a[0][3] = -det3(a30) / d;
  a[1][0] = -det3(a01) / d;
  a[1][1] = det3(a11) / d;
  a[1][2] = -det3(a21) / d;
  a[1][3] = det3(a31) / d;
  a[2][0] = det3(a02) / d;
  a[2][1] = -det3(a12) / d;
  a[2][2] = det3(a22) / d;
  a[2][3] = -det3(a32) / d;
  a[3][0] = -det3(a03) / d;
  a[3][1] = det3(a13) / d;
  a[3][2] = -det3(a23) / d;
  a[3][3] = det3(a33) / d;

  return a;
}
function inverse(m) {
  if (m.matrix != true) console.log("not a matrix");
  if (m.length == 2) return inverse2(m);
  if (m.length == 3) return inverse3(m);
  if (m.length == 4) return inverse4(m);
}

function normalMatrix(m, flag) {
  var a = mat4();
  a = inverse(transpose(m));
  if (flag != true) return a;else {
    var b = mat3();
    for (var i = 0; i < 3; i++) {
      for (var j = 0; j < 3; j++) {
        b[i][j] = a[i][j];
      }
    }return b;
  }
}

;!function ($) {
  $.fn.classes = function (callback) {
    var classes = [];
    $.each(this, function (i, v) {
      var splitClassName = v.className.split(/\s+/);
      for (var j = 0; j < splitClassName.length; j++) {
        var className = splitClassName[j];
        if (-1 === classes.indexOf(className)) {
          classes.push(className);
        }
      }
    });
    if ('function' === typeof callback) {
      for (var i in classes) {
        callback(classes[i]);
      }
    }
    return classes;
  };
}(jQuery);

var initGlProgram = function initGlProgram(gl, shaders) {
  var program = gl.createProgram();
  for (var i = 0; i < shaders.length; i += 1) {
    gl.attachShader(program, shaders[i]);
  }
  gl.linkProgram(program);

  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    window.console.error("failed to link: " + gl.getProgramInfoLog(program));
    return null;
  }
  return program;
};

var initGlShader = function initGlShader(gl, shaderScript, shaderType) {
  var shader = gl.createShader(shaderType);
  gl.shaderSource(shader, shaderScript);
  gl.compileShader(shader);
  return shader;
};

var BgSprite = function () {
  function BgSprite(editor) {
    _classCallCheck(this, BgSprite);

    this.map = editor.map;
    this.style = editor.style;
    this.editor = editor;
    this.animating = false;
    this.options = {
      timeGap: 100,
      duration: 1000,
      entryLength: 100
    };
    this.width = this.editor.canvasNode.width;
    this.height = this.editor.canvasNode.height;
    this.ctx = this.editor.bgCtx;
    this.gl = this.editor.bgGl;
    this.init();
  }

  _createClass(BgSprite, [{
    key: "init",
    value: function init() {}
  }, {
    key: "update",
    value: function update() {
      this.map = this.editor.map;
      this.style = this.editor.style;
      this.width = this.editor.canvasNode.width;
      this.height = this.editor.canvasNode.height;
      this.ctx = this.editor.bgCtx;
      this.gl = this.editor.bgGl;
    }
  }, {
    key: "drawStatic",
    value: function drawStatic() {}
  }, {
    key: "drawFrame",
    value: function drawFrame() {}
  }, {
    key: "advance",
    value: function advance(t) {
      var stop = true;
      return stop;
    }
  }, {
    key: "clear",
    value: function clear() {}
  }]);

  return BgSprite;
}();

var TextSprite = function () {
  function TextSprite(editor, options) {
    _classCallCheck(this, TextSprite);

    this.map = editor.map;
    this.style = editor.style;
    this.editor = editor;
    this.animating = false;
    this.textStyle = 1;
    this.options = options || {
      timeGap: 100,
      duration: 2000,
      entryLength: 100
    };
    this.ctx = this.editor.ctx;
    this.gl = this.editor.gl;
    this.animationStyle = 0;
    this._t = 0.0;
    this.init();
  }

  _createClass(TextSprite, [{
    key: "init",
    value: function init() {
      this.clearArea = [0, 0, this.editor.ctx.canvas.width, this.editor.ctx.canvas.height];
    }
  }, {
    key: "changeOptions",
    value: function changeOptions(options) {
      Object.keys(options).forEach(function (key) {
        this.options[key] = options[key];
      });
    }
  }, {
    key: "changeTextStyle",
    value: function changeTextStyle(style) {
      this.textStyle = style;
    }
  }, {
    key: "update",
    value: function update() {
      this.map = this.editor.map;
      this.style = this.editor.style;
      this.startTime = Date.now();
    }
  }, {
    key: "drawStatic",
    value: function drawStatic() {}
  }, {
    key: "drawFrame",
    value: function drawFrame() {}
  }, {
    key: "advance",
    value: function advance(t) {
      var stop = true;
      if (stop) {
        this.stop();
      }
      return stop;
    }
  }, {
    key: "stop",
    value: function stop() {
      this._t = 0;
      this.animating = false;
    }
  }, {
    key: "getAnimationMap",
    value: function getAnimationMap() {
      if (!this.animationMap) {
        this.animationMap = [{
          label: '无',
          value: 0
        }];
      }
      return this.animationMap;
    }
  }, {
    key: "clear",
    value: function clear(useScissor) {
      if (useScissor) {
        if (!this.clearArea.length) {
          return;
        }
        this.editor.ctx.clearRect(this.clearArea[0], this.clearArea[1], this.clearArea[2], this.clearArea[3]);
      }
      this.editor.ctx.clearRect(0, 0, this.editor.ctx.canvas.width, this.editor.ctx.canvas.height);
    }
  }]);

  return TextSprite;
}();

var util = {
  hexToRgb: function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
      r: parseInt(result[1], 16),
      g: parseInt(result[2], 16),
      b: parseInt(result[3], 16)
    } : null;
  },

  rgbToHex: function rgbToHex(r, g, b) {
    function componentToHex(c) {
      var hex = c.toString(16);
      return hex.length == 1 ? "0" + hex : hex;
    }

    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
  },

  getSpriteEntity: function () {
    var entities = [];
    return function (className, editor) {
      var Klass = eval(className);
      return entities[className] ? entities[className] : entities[className] = new Klass(editor);
    };
  }(),

  getInitialState: function getInitialState(name) {
    return config[name + 'Map'][config.state[name + 'Index']].value;
  }
};
var config = {
  fontMap: [{
    label: '隶书',
    value: '隶书'
  }, {
    label: '楷体',
    value: '楷体'
  }, {
    label: '华文新魏',
    value: '华文新魏'
  }],

  fontSizeMap: [{
    label: '小',
    value: 40
  }, {
    label: '中',
    value: 50
  }, {
    label: '大',
    value: 60
  }],

  fontColorMap: [{
    label: '黑',
    value: '#000000'
  }, {
    label: '白',
    value: '#FFFFFF'
  }, {
    label: '褐',
    value: '#514848'
  }],

  textStyleMap: [{
    Klass: 'GlTextSprite',
    style: 1,
    label: '渲墨',
    value: 0
  }, {
    Klass: 'GlTextSprite',
    style: 2,
    label: '立体',
    value: 1
  }],

  textAlignMap: [{
    label: '居中',
    value: 'center'
  }, {
    label: '左对齐',
    value: 'left'
  }],

  backgroundMap: [{
    Klass: 'PureBgSprite',
    label: '纯色',
    value: 0,
    colors: ['rgb(235, 235, 235)', '#FEFEFE', '#3a3a3a']
  }, {
    Klass: 'TreeBgSprite',
    label: '月下林间',
    value: 1,
    colors: ['rgb(235, 235, 235)', '#b1a69b', '#3a3a3a']
  }],

  animationMap: [{
    label: '无',
    value: 0
  }, {
    label: '粒子缓入',
    value: 1
  }],

  state: {
    fontIndex: 0,
    fontSizeIndex: 0,
    fontColorIndex: 0,
    textStyleIndex: 0,
    textAlignIndex: 0,
    backgroundIndex: 0,
    animationIndex: 1,
    bgColorIndex: 0
  }
};

var textCtx = document.createElement("canvas").getContext("2d");

function calTriangleArea(A, B, C) {
  var a = distance(A, B);
  var b = distance(B, C);
  var c = distance(C, A);

  var p = (a + b + c) / 2;
  return Math.sqrt(p * (p - a) + (p - b) + (p - c));
}

function distance(p1, p2) {
  return Math.sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
}

function makeTextData(text, size, options) {
  var pointsArray = [];
  var normalsArray = [];
  var colorsArray = [];
  var compIdxArray = [];
  var e = [];
  var g = [];
  var markMap = [];
  var cnt = 0;
  var width = size;
  var height = size;
  textCtx.canvas.width = width;
  textCtx.canvas.height = height;

  var ctx = textCtx;
  ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
  ctx.font = size + 'px ' + (options.font || '隶书');
  ctx.fillStyle = 'black';
  ctx.strokeStyle = 'black';
  ctx.lineWidth = 1;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(text, width / 2, height / 2);
  var imageData = ctx.getImageData(0, 0, ctx.canvas.width, ctx.canvas.height);
  var data = imageData.data;
  console.log(data);
  var grid = 1;

  var triangle = function triangle(a, b, c, index, options) {
    pointsArray.push(a);
    pointsArray.push(b);
    pointsArray.push(c);

    compIdxArray.push(index);
    compIdxArray.push(index);
    compIdxArray.push(index);

    var normal = cross(subtract(b, a), subtract(c, a));
    normalsArray.push(vec4(normal[0], normal[1], normal[2], 0.0));
    normalsArray.push(vec4(normal[0], normal[1], normal[2], 0.0));
    normalsArray.push(vec4(normal[0], normal[1], normal[2], 0.0));

    if (options) {
      switch (options.colorStyle) {
        case 'test':
          var index = options.index;
          var n = options.n;
          var area = calTriangleArea(a, b, c);

          var color = vec4(0.0, 0.0, 0.0, area > 100 ? 1.0 : Math.pow(area / 100, 0.5));
          colorsArray.push(color);
          colorsArray.push(color);
          colorsArray.push(color);
          break;
      }
    }
  };

  var hasPixel = function hasPixel(j, i) {
    //第j行，第i列
    if (i < 0 || j < 0) {
      return false;
    }
    return !!data[(j * ctx.canvas.width + i) * 4 + 3];
  };

  var markPoint = function markPoint(j, i) {
    var value = 0;

    if (i > 0 && hasPixel(j, i - 1)) {
      //与左边连通
      value = g[j][i - 1];
    } else {
      value = ++cnt;
    }

    if (j > 0 && hasPixel(j - 1, i) && (i === 0 || !hasPixel(j - 1, i - 1))) {
      //与上连通 且 与左上不连通 （即首次和上一行连接）
      if (g[j - 1][i] !== value) {
        markMap.push([g[j - 1][i], value]);
      }
    }

    if (!hasPixel(j, i - 1)) {
      //行首
      if (hasPixel(j - 1, i - 1) && g[j - 1][i - 1] !== value) {
        //与左上连通
        markMap.push([g[j - 1][i - 1], value]);
      }
    }

    if (!hasPixel(j, i + 1)) {
      //行尾
      if (hasPixel(j - 1, i + 1) && g[j - 1][i + 1] !== value) {
        //与右上连通
        markMap.push([g[j - 1][i + 1], value]);
      }
    }

    return value;
  };

  for (var j = 0; j < ctx.canvas.height; j += grid) {
    g.push([]);
    e.push([]);

    for (var i = 0; i < ctx.canvas.width; i += grid) {
      var value = 0;
      var isEdge = false;

      if (hasPixel(j, i)) {
        value = markPoint(j, i);
      }
      e[j][i] = isEdge;
      g[j][i] = value;
    }
  }

  var finalGraph = seperateGraph(g, e, markMap, cnt);
  var graphs = finalGraph[0];
  var sampledEdge = finalGraph[2];
  var graphMap = [];
  var z = options.thick / 2;
  var z2 = -z;

  sampledEdge.forEach(function (graph, index) {
    var points = graph[0];
    var holes = graph[1];
    var triangles = earcut(flatten(points), holes);
    var num = points.length;

    for (var n = 0; n < triangles.length; n += 3) {
      var a = points[triangles[n]];
      var b = points[triangles[n + 1]];
      var c = points[triangles[n + 2]];

      //=====字体正面数据=====
      triangle(vec3(a[0], a[1], z), vec3(b[0], b[1], z), vec3(c[0], c[1], z), index);

      //=====字体背面数据=====
      triangle(vec3(a[0], a[1], z2), vec3(b[0], b[1], z2), vec3(c[0], c[1], z2), index);
    }

    var holesMap = [];
    var last = 0;

    if (holes.length) {
      for (var holeIndex = 0; holeIndex < holes.length; holeIndex++) {
        holesMap.push([last, holes[holeIndex] - 1]);
        last = holes[holeIndex];
      }
    }

    holesMap.push([last, points.length - 1]);

    for (var i = 0; i < holesMap.length; i++) {
      var startAt = holesMap[i][0];
      var endAt = holesMap[i][1];

      for (var j = startAt; j < endAt; j++) {
        triangle(vec3(points[j][0], points[j][1], z), vec3(points[j][0], points[j][1], z2), vec3(points[j + 1][0], points[j + 1][1], z), index);
        triangle(vec3(points[j][0], points[j][1], z2), vec3(points[j + 1][0], points[j + 1][1], z2), vec3(points[j + 1][0], points[j + 1][1], z), index);
      }
      triangle(vec3(points[startAt][0], points[startAt][1], z), vec3(points[startAt][0], points[startAt][1], z2), vec3(points[endAt][0], points[endAt][1], z), index);
      triangle(vec3(points[startAt][0], points[startAt][1], z2), vec3(points[endAt][0], points[endAt][1], z2), vec3(points[endAt][0], points[endAt][1], z), index);
    }
    graphMap.push(pointsArray.length);
  });
  return [pointsArray, normalsArray, colorsArray, compIdxArray, graphMap];
}

function seperateGraph(g, e, markMap, cnt) {
  var graphs = [];
  var graphsEdge = [];
  var sampledEdge = [];

  var markArr = [];
  var fathers = [];
  var markCollection = [];

  for (var i = 0; i < cnt; i++) {
    markArr[i] = i;
  }

  var findFather = function findFather(n) {
    if (markArr[n] === n) {
      return n;
    } else {
      markArr[n] = findFather(markArr[n]);
      return markArr[n];
    }
  };

  for (i = 0; i < markMap.length; i++) {
    var a = markMap[i][0];
    var b = markMap[i][1];

    var f1 = findFather(a);
    var f2 = findFather(b);

    if (f1 !== f2) {
      markArr[f2] = f1;
    }
  }

  for (i = 1; i < markArr.length; i++) {
    var f = findFather(markArr[i]);
    var index = fathers.indexOf(f);
    if (index !== -1) {
      markCollection[index].push(i);
    } else {
      fathers.push(f);
      markCollection.push([i]);
    }
  }

  for (i = 0; i < markCollection.length; i++) {
    graphs.push([]);
    graphsEdge.push([]);
  }

  for (var j = 0; j < g.length; j++) {
    for (i = 0; i < g[j].length; i++) {
      var v = g[j][i];
      for (var n = 0; n < markCollection.length; n++) {
        if (markCollection[n].indexOf(v) !== -1) {
          g[j][i] = n + 1;
          var p = [i, j];
          graphs[n].push(p);
          if (e[j][i]) {
            graphsEdge[n].push(p);
          }
        }
      }
    }
  }

  graphs.forEach(function (shape, index) {
    sampledEdge.push(orderEdge(g, e, index, 0));
  });
  return [graphs, graphsEdge, sampledEdge];
}

function findOuterContourEntry(g, v) {
  var start = [-1, -1];
  for (var j = 0; j < g.length; j++) {
    for (var i = 0; i < g[0].length; i++) {
      if (g[j][i] === v) {
        start = [i, j];
        return start;
      }
    }
  }
  return start;
}

function findInnerContourEntry(g, v, e) {
  var start = false;
  for (var j = 0; j < g.length; j++) {
    for (var i = 0; i < g[0].length; i++) {
      if (g[j][i] === v && g[j + 1] && g[j + 1][i] === 0) {
        var isInContours = false;
        if (typeof e[j][i] === 'number') {
          isInContours = true;
        }
        if (!isInContours) {
          start = [i, j];
          return start;
        }
      }
    }
  }
  return start;
}

function orderEdge(g, e, v, gap) {
  v++;
  var rs = [];
  var entryRecord = [];
  var start = findOuterContourEntry(g, v);
  var next = start;
  var end = false;
  rs.push(start);
  entryRecord.push(6);
  var holes = [];
  var mark;
  var holeMark = 2;
  e[start[1]][start[0]] = holeMark;

  var process = function process(i, j) {
    if (i < 0 || i >= g[0].length || j < 0 || j >= g.length) {
      return false;
    }

    if (g[j][i] !== v || tmp) {
      return false;
    }

    e[j][i] = holeMark;
    tmp = [i, j];
    rs.push(tmp);
    mark = true;

    return true;
  };

  var map = [function (i, j) {
    return { 'i': i + 1, 'j': j };
  }, function (i, j) {
    return { 'i': i + 1, 'j': j + 1 };
  }, function (i, j) {
    return { 'i': i, 'j': j + 1 };
  }, function (i, j) {
    return { 'i': i - 1, 'j': j + 1 };
  }, function (i, j) {
    return { 'i': i - 1, 'j': j };
  }, function (i, j) {
    return { 'i': i - 1, 'j': j - 1 };
  }, function (i, j) {
    return { 'i': i, 'j': j - 1 };
  }, function (i, j) {
    return { 'i': i + 1, 'j': j - 1 };
  }];

  var convertEntry = function convertEntry(index) {
    var arr = [4, 5, 6, 7, 0, 1, 2, 3];
    return arr[index];
  };

  while (!end) {
    var i = next[0];
    var j = next[1];
    var tmp = null;
    var entryIndex = entryRecord[entryRecord.length - 1];

    for (var c = 0; c < 8; c++) {
      var index = (entryIndex + 1 + c) % 8;
      var hasNext = process(map[index](i, j).i, map[index](i, j).j);
      if (hasNext) {
        entryIndex = convertEntry(index);
        break;
      }
    }

    if (tmp) {
      next = tmp;

      if (next[0] === start[0] && next[1] === start[1]) {
        var innerEntry = findInnerContourEntry(g, v, e);
        if (innerEntry) {
          next = start = innerEntry;
          e[start[1]][start[0]] = holeMark;
          rs.push(next);
          entryRecord.push(entryIndex);
          entryIndex = 2;
          holes.push(rs.length - 1);
          holeMark++;
        } else {
          end = true;
        }
      }
    } else {
      rs.splice(rs.length - 1, 1);
      entryIndex = convertEntry(entryRecord.splice(entryRecord.length - 1, 1)[0]);
      next = rs[rs.length - 1];
    }

    entryRecord.push(entryIndex);
  }
  return [rs, holes];
}

function convertCanvasToGl(x, y, z, size) {
  return vec3((x - size / 2) / size, (size / 2 - y) / size, z);
}

function convertCanvasToGl2D(x, y, size) {
  return vec2((x - size / 2) / size, (size / 2 - y) / size);
}

function parseArr(data, size) {
  var points = data[0];
  var normals = data[1];
  var colors = data[2];

  for (var i = 0; i < points.length; i++) {
    points[i] = convertCanvasToGl(points[i][0], points[i][1], points[i][2], size);
  }
  return [points, normals, colors, data[3], data[4]];
}

function hslToRgb(h, s, l) {
  var r, g, b;

  if (s == 0) {
    r = g = b = l;
  } else {
    var hue2rgb = function hue2rgb(p, q, t) {
      if (t < 0) t += 1;
      if (t > 1) t -= 1;
      if (t < 1 / 6) return p + (q - p) * 6 * t;
      if (t < 1 / 2) return q;
      if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
      return p;
    };

    var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    var p = 2 * l - q;
    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }

  return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
}

function rgbToHsl(r, g, b) {
  r /= 255, g /= 255, b /= 255;
  var max = Math.max(r, g, b),
      min = Math.min(r, g, b);
  var h,
      s,
      l = (max + min) / 2;

  if (max == min) {
    h = s = 0;
  } else {
    var d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);break;
      case g:
        h = (b - r) / d + 2;break;
      case b:
        h = (r - g) / d + 4;break;
    }
    h /= 6;
  }

  return [h, s, l];
}

function componentToHex(c) {
  var hex = c.toString(16);
  return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
  return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

function hexToRgb(hex) {
  var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result ? {
    r: parseInt(result[1], 16),
    g: parseInt(result[2], 16),
    b: parseInt(result[3], 16)
  } : null;
}

var GlTextSprite = function (_TextSprite) {
  _inherits(GlTextSprite, _TextSprite);

  function GlTextSprite(editor, options) {
    _classCallCheck(this, GlTextSprite);

    return _possibleConstructorReturn(this, (GlTextSprite.__proto__ || Object.getPrototypeOf(GlTextSprite)).call(this, editor, options));
  }

  _createClass(GlTextSprite, [{
    key: "init",
    value: function init() {
      _get(GlTextSprite.prototype.__proto__ || Object.getPrototypeOf(GlTextSprite.prototype), "init", this).call(this);
      this.canvas = this.editor.glcanvasNode;
      this.dataMap = {};
      this.initGl();
    }
  }, {
    key: "update",
    value: function update() {
      _get(GlTextSprite.prototype.__proto__ || Object.getPrototypeOf(GlTextSprite.prototype), "update", this).call(this);
      this._finish = false;
      this.clearArea = [0, 0, this.canvas.width, this.canvas.height];
    }
  }, {
    key: "parseData",
    value: function parseData() {
      var map = this.map;
      for (var i = 0; i < map.length; i++) {
        for (var j = 0; j < map[i].length; j++) {
          if (map[i][j].char && map[i][j].char !== ' ') {
            var data = this.getTextData(map[i][j].char);
            map[i][j].data = data;
          }
        }
      }
      $('.render-tip').removeClass('show');
    }
  }, {
    key: "getTextData",
    value: function getTextData(char) {
      var font = this.style.font;
      var data;

      if (!this.dataMap[font]) {
        this.dataMap[font] = {};
      }

      if (Object.keys(this.dataMap[font]).indexOf(char) === -1) {
        var size = 200;
        data = makeTextData(char, size, {
          thick: 0.1,
          font: font
        });
        data = parseArr(data, size);
        this.dataMap[font][char] = data;
      } else {
        data = this.dataMap[font][char];
      }
      return data;
    }
  }, {
    key: "drawStatic",
    value: function drawStatic() {
      var ctx = this.editor.ctx;
      var map = this.map;

      this.parseData();
      this.clear();

      for (var i = 0; i < map.length; i++) {
        for (var j = 0; j < map[i].length; j++) {
          if (map[i][j].char && map[i][j].char !== ' ') {
            var data = map[i][j].data;
            this.renderGlText(data, {
              pos: {
                x: map[i][j].x + this.editor.style.fontSize / 2,
                y: map[i][j].y + this.style.lineHeight() / 2
              },
              fontSize: this.style.fontSize,
              fontColor: this.style.fontColor,
              textStyle: this.textStyle,
              t: 1.0
            });
          }
        }
      }
    }
  }, {
    key: "drawFrame",
    value: function drawFrame() {
      switch (this.style.animation) {
        case 0:
          var map = this.map;
          var fontSize = this.editor.style.fontSize;
          var gl = this.gl;
          var finish = true;

          this.clear(true);

          for (var i = 0; i < map.length; i++) {
            for (var j = 0; j < map[i].length; j++) {
              if (map[i][j].char && map[i][j].char !== ' ') {
                var data = map[i][j].data;
                var period = this.t;
                var t = 0;

                if (Math.floor(period / 1000) < i) {
                  map[i][j].finish = -1;
                  finish = false;
                } else if (Math.floor(period / 1000) === i) {
                  map[i][j].finish = 0;
                  finish = false;
                  t = period % 1000 / 1000;
                } else if (Math.floor(period / 1000) > i) {
                  t = 1.0;
                  if (map[i][j].finish) {
                    map[i][j].finish += 1;
                    this.clearArea = [0, 0, this.canvas.width, this.canvas.height - (map[i][j].y + this.style.lineHeight())];
                  } else {
                    finish = false;
                    map[i][j].finish = 1;
                  }
                }

                if (map[i][j].finish >= 0 && map[i][j].finish <= 2) {
                  var pos = {
                    x: map[i][j].x + this.editor.style.fontSize / 2,
                    y: map[i][j].y + this.style.lineHeight() / 2
                  };
                  this.renderGlText(data, {
                    pos: pos,
                    fontSize: fontSize,
                    fontColor: this.style.fontColor,
                    textStyle: this.textStyle,
                    t: t
                  });
                }
              }
            }
          }

          this._finish = finish;
          break;
        case 1:
          this.drawStatic();
          break;
      }
    }
  }, {
    key: "advance",
    value: function advance(t) {
      var stop = false;
      this.t = t;

      switch (this.style.animation) {
        case 0:
          stop = this._finish;
          break;
        case 1:
          stop = true;
          break;
      }

      if (stop) {
        this.stop();
      }

      return stop;
    }
  }, {
    key: "stop",
    value: function stop() {
      this._t = 0;
      this.animating = false;
      this.clearArea = [0, 0, this.canvas.width, this.canvas.height];
    }
  }, {
    key: "getAnimationMap",
    value: function getAnimationMap() {
      if (!this.animationMap) {
        this.animationMap = [{
          label: '纵横捭阖',
          value: 0
        }, {
          label: '无',
          value: 1
        }];
      }
      return this.animationMap;
    }
  }, {
    key: "initGl",
    value: function initGl() {
      var gl = this.gl;

      this.vertexShader = initGlShader(gl, "attribute vec3 vPosition;\n      attribute vec4 vNormal;\n      attribute float compIdx;\n\n      varying vec4 fColor;\n\n      uniform float t;\n      uniform float compNum;\n      uniform int style;\n      uniform vec3 vColor;\n\n      uniform vec4 ambientProduct, diffuseProduct, specularProduct;\n      uniform mat4 modelViewMatrix;\n      uniform mat4 projectionMatrix;\n      uniform vec4 lightPosition;\n      uniform float shininess;\n      uniform mat3 normalMatrix;\n\n      void main() {\n        float n = (1.0 - t) * pow((1.0 + compIdx) / compNum, 2.0);\n        float xx = vPosition.x - vPosition.x * pow(n, 0.5);\n        float yy = vPosition.y - 2.0 * vPosition.y * pow(n, 2.0);\n        float zz = vPosition.z - pow(n, 0.5);\n        float ww = 1.0;\n        \n        vec4 aPosition = vec4(xx, yy, zz, ww);\n        vec3 pos = (modelViewMatrix * aPosition).xyz;\n\n        vec3 L;\n        \n        if(lightPosition.w == 0.0) \n          L = normalize(lightPosition.xyz);\n        else \n          L = normalize( lightPosition.xyz - pos );\n\n        vec3 E = -normalize( pos );\n        vec3 H = normalize( L + E );\n        vec3 N = normalize( normalMatrix*vNormal.xyz);\n\n        vec4 ambient = vec4(vColor, 1.0);\n\n        float Kd = max( dot(L, N), 0.0 );\n        vec4  diffuse = Kd*diffuseProduct;\n\n        float Ks = pow( max(dot(N, H), 0.0), shininess );\n        vec4  specular = Ks * specularProduct;\n        \n        if( dot(L, N) < 0.0 ) {\n          specular = vec4(0.0, 0.0, 0.0, 1.0);\n        }\n\n        gl_PointSize = 1.0;\n        gl_Position = projectionMatrix * modelViewMatrix * aPosition;\n\n        if (style == 1) {\n          if (dot(N, L) < 0.0) {\n            fColor = vec4(vColor, 1.0);\n          } else if (dot(N, L) < 0.8) {\n            fColor = vec4(0.0, 0.0, 0.0, 0.0);\n          } else {\n            fColor = vec4(vColor, 1.0);\n          }\n        } else if (style == 2) {\n          fColor = ambient + diffuse +specular;\n        }\n      }", gl.VERTEX_SHADER);

      this.fragmentShader = initGlShader(gl, "precision mediump float;\n\n      varying vec4 fColor;\n\n      void main() {\n        gl_FragColor = fColor;\n      }", gl.FRAGMENT_SHADER);

      this.program = initGlProgram(gl, [this.vertexShader, this.fragmentShader]);
    }
  }, {
    key: "renderGlText",
    value: function renderGlText(data, options) {
      var pos = options.pos;
      var fontSize = options.fontSize;
      var fontColor = hexToRgb(options.fontColor);
      var textStyle = options.textStyle;
      fontColor = vec3(fontColor.r / 255, fontColor.g / 255, fontColor.b / 255);

      var t = options.t;

      if (t === 0) {
        return;
      }

      var gl = this.gl;
      var program = this.program;

      var pointsArray = [];
      var normalsArray = [];
      var colorArray = [];
      var compIdxArray = [];

      var near = -10;
      var far = 10;
      var radius = 1.5;
      var theta = 0.7;
      var phi = 0.1;
      var dr = 5.0 * Math.PI / 180.0;

      var left = -3.0;
      var right = 3.0;
      var ytop = 3.0;
      var bottom = -3.0;

      var lightPosition = vec4(1.0, 1.0, 1.0, 0.0);
      var lightAmbient = vec4(0.5, 0.5, 0.5, 1.0);
      var lightDiffuse = vec4(1.0, 1.0, 1.0, 1.0);
      var lightSpecular = vec4(1.0, 1.0, 1.0, 1.0);

      var materialAmbient = vec4(0.5, 0.5, 0.5, 1.0);
      var materialDiffuse = vec4(0.7, 0.7, 0.7, 1.0);
      var materialSpecular = vec4(1.0, 1.0, 1.0, 1.0);
      var materialShininess = 20.0;

      var ambientColor, diffuseColor, specularColor;

      var modelViewMatrix, projectionMatrix;
      var modelViewMatrixLoc, projectionMatrixLoc;

      var normalMatrix, normalMatrixLoc;

      var eye;
      var at = vec3(0.0, 0.0, 0.0);
      var up = vec3(0.0, 1.0, 0.0);

      gl.useProgram(program);
      gl.enable(gl.DEPTH_TEST);

      gl.viewport(pos.x - fontSize * 4, this.canvas.height - pos.y - fontSize * 4, fontSize * 8, fontSize * 8);

      var ambientProduct = mult(lightAmbient, materialAmbient);
      var diffuseProduct = mult(lightDiffuse, materialDiffuse);
      var specularProduct = mult(lightSpecular, materialSpecular);
      var data;

      pointsArray = data[0];
      normalsArray = data[1];
      colorArray = data[2];
      var compIdxArray = data[3];
      var graphMap = data[4];

      var nBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, nBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(normalsArray), gl.STATIC_DRAW);

      var vNormal = gl.getAttribLocation(program, "vNormal");
      gl.vertexAttribPointer(vNormal, 4, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(vNormal);

      var vBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(pointsArray), gl.STATIC_DRAW);

      var vPosition = gl.getAttribLocation(program, "vPosition");
      gl.vertexAttribPointer(vPosition, 3, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(vPosition);

      var iBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, iBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(compIdxArray), gl.STATIC_DRAW);

      var compIdx = gl.getAttribLocation(program, "compIdx");
      gl.vertexAttribPointer(compIdx, 1, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(compIdx);

      modelViewMatrixLoc = gl.getUniformLocation(program, "modelViewMatrix");
      projectionMatrixLoc = gl.getUniformLocation(program, "projectionMatrix");
      normalMatrixLoc = gl.getUniformLocation(program, "normalMatrix");

      gl.uniform4fv(gl.getUniformLocation(program, "ambientProduct"), flatten(ambientProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "diffuseProduct"), flatten(diffuseProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "specularProduct"), flatten(specularProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "lightPosition"), flatten(lightPosition));
      gl.uniform1f(gl.getUniformLocation(program, "shininess"), materialShininess);
      gl.uniform1f(gl.getUniformLocation(program, "compNum"), graphMap.length);
      gl.uniform1f(gl.getUniformLocation(program, "t"), t);
      gl.uniform1i(gl.getUniformLocation(program, "style"), textStyle);
      gl.uniform3fv(gl.getUniformLocation(program, "vColor"), flatten(fontColor));

      eye = vec3(radius * Math.sin(theta) * Math.cos(phi), radius * Math.sin(theta) * Math.sin(phi), radius * Math.cos(theta));

      modelViewMatrix = lookAt(eye, at, up);
      projectionMatrix = ortho(left, right, bottom, ytop, near, far);

      normalMatrix = [vec3(modelViewMatrix[0][0], modelViewMatrix[0][1], modelViewMatrix[0][2]), vec3(modelViewMatrix[1][0], modelViewMatrix[1][1], modelViewMatrix[1][2]), vec3(modelViewMatrix[2][0], modelViewMatrix[2][1], modelViewMatrix[2][2])];
      gl.uniformMatrix4fv(modelViewMatrixLoc, false, flatten(modelViewMatrix));
      gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix));
      gl.uniformMatrix3fv(normalMatrixLoc, false, flatten(normalMatrix));
      gl.drawArrays(gl.TRIANGLES, 0, pointsArray.length);
    }
  }, {
    key: "clear",
    value: function clear(useScissor) {
      var gl = this.gl;
      if (useScissor) {
        if (!this.clearArea.length) {
          return;
        }
        gl.enable(gl.SCISSOR_TEST);
        gl.scissor(this.clearArea[0], this.clearArea[1], this.clearArea[2], this.clearArea[3]);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.disable(gl.SCISSOR_TEST);
      } else {
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      }
    }
  }]);

  return GlTextSprite;
}(TextSprite);

var PureBgSprite = function (_BgSprite) {
  _inherits(PureBgSprite, _BgSprite);

  function PureBgSprite(editor) {
    _classCallCheck(this, PureBgSprite);

    return _possibleConstructorReturn(this, (PureBgSprite.__proto__ || Object.getPrototypeOf(PureBgSprite)).call(this, editor));
  }

  _createClass(PureBgSprite, [{
    key: "drawStatic",
    value: function drawStatic() {
      this.clear();
      var ctx = this.ctx;
      var color = config.backgroundMap[config.state.backgroundIndex].colors[config.state.bgColorIndex];
      ctx.fillStyle = color;
      ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    }
  }, {
    key: "drawFrame",
    value: function drawFrame() {
      this.drawStatic();
    }
  }, {
    key: "advance",
    value: function advance(t) {
      var stop = true;
      return stop;
    }
  }, {
    key: "clear",
    value: function clear() {
      this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
      this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
    }
  }]);

  return PureBgSprite;
}(BgSprite);

var TreeBgSprite = function (_BgSprite2) {
  _inherits(TreeBgSprite, _BgSprite2);

  function TreeBgSprite(editor) {
    _classCallCheck(this, TreeBgSprite);

    return _possibleConstructorReturn(this, (TreeBgSprite.__proto__ || Object.getPrototypeOf(TreeBgSprite)).call(this, editor));
  }

  _createClass(TreeBgSprite, [{
    key: "init",
    value: function init() {
      this.initGl();
    }
  }, {
    key: "drawStatic",
    value: function drawStatic() {
      this.clear();
      var color = config.backgroundMap[config.state.backgroundIndex].colors[config.state.bgColorIndex];
      this.ctx.fillStyle = color;
      this.ctx.fillRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);

      this.renderGlMoon({
        pos: { x: 150, y: 150 },
        r: 100,
        t: 1.0,
        numTimesToSubdivide: 5
      });

      this.renderGlTree({
        pos: { x: this.ctx.canvas.width - 350, y: this.ctx.canvas.height },
        theta: 0.2,
        minH: this.ctx.canvas.width / 80,
        iniH: 30
      });

      this.renderGlTree({
        pos: { x: this.ctx.canvas.width - 200, y: this.ctx.canvas.height },
        theta: 0.2,
        minH: this.ctx.canvas.width / 65,
        iniH: 50
      });
    }
  }, {
    key: "drawFrame",
    value: function drawFrame() {
      this.clear();
      var color = config.backgroundMap[config.state.backgroundIndex].colors[config.state.bgColorIndex];
      this.ctx.fillStyle = color;
      this.ctx.fillRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);

      this.renderGlMoon({
        pos: { x: 150, y: 150 },
        r: 100,
        t: this._t,
        numTimesToSubdivide: 5
      });

      this.renderGlTree({
        pos: { x: this.ctx.canvas.width - 350, y: this.ctx.canvas.height },
        theta: 0.2 * this._t,
        minH: this.ctx.canvas.width / 80,
        iniH: 30
      });

      this.renderGlTree({
        pos: { x: this.ctx.canvas.width - 200, y: this.ctx.canvas.height },
        theta: 0.2 * this._t,
        minH: this.ctx.canvas.width / 65,
        iniH: 50
      });
    }
  }, {
    key: "advance",
    value: function advance(t) {
      var stop = true;
      this.t = t;

      if (this._t < 1) {
        var dt = this.options.duration;
        var progress = t / dt;
        this._t = Math.pow(progress, 1);
        stop = false;
      }

      if (stop) {
        this.stop();
      }

      return stop;
    }
  }, {
    key: "stop",
    value: function stop() {
      this.animating = false;
      this._t = 0;
    }
  }, {
    key: "initGl",
    value: function initGl() {
      var gl = this.gl;

      this.vertexMoonShader = initGlShader(gl, "attribute vec4 vPosition;\n      attribute vec4 vNormal;\n      attribute float compIdx;\n\n      varying vec4 fColor;\n\n      uniform float t;\n      uniform float compNum;\n      uniform vec3 vColor;\n\n      uniform vec4 ambientProduct, diffuseProduct, specularProduct;\n      uniform mat4 modelViewMatrix;\n      uniform mat4 projectionMatrix;\n      uniform vec4 lightPosition;\n      uniform float shininess;\n      uniform mat3 normalMatrix;\n\n      void main() {\n        vec3 pos = (modelViewMatrix * vPosition).xyz;\n\n        vec3 L;\n        \n        if(lightPosition.w == 0.0) L = normalize(lightPosition.xyz);\n        else L = normalize( lightPosition.xyz - pos );\n         \n        vec3 E = -normalize( pos );\n        \n        vec3 H = normalize( L + E );\n        \n        vec3 N = normalize( normalMatrix*vNormal.xyz);\n\n        vec4 ambient = ambientProduct;\n\n        float Kd = max( dot(L, N), 0.0 );\n        vec4  diffuse = Kd*diffuseProduct;\n\n        float Ks = pow( max(dot(N, H), 0.0), shininess );\n        vec4  specular = Ks * specularProduct;\n        \n        if( dot(L, N) < 0.0 ) {\n          specular = vec4(0.0, 0.0, 0.0, 1.0);\n        }\n\n        gl_PointSize = 1.0;\n        gl_Position = projectionMatrix * modelViewMatrix * vPosition;\n\n        fColor = ambient + diffuse +specular;\n        fColor.a = 1.0;\n\n        if (dot(N, L) < 0.0) {\n          fColor = vec4(1.0, 1.0, 0.9 - 0.2 * t, 1.0);\n        } else if (dot(N, L) < 0.9 + t * 0.09 ) {\n          fColor = vec4(1.0, 1.0, 0.9 - 0.1 * t, 1.0);\n        } else {\n          fColor = vec4(1.0, 1.0, 0.9 + 0.05 * t, 1.0);\n        }\n      }", gl.VERTEX_SHADER);

      this.fragmentMoonShader = initGlShader(gl, "precision mediump float;\n\n      varying vec4 fColor;\n\n      void main() {\n        gl_FragColor = fColor;\n      }", gl.FRAGMENT_SHADER);

      this.moonProgram = initGlProgram(gl, [this.vertexMoonShader, this.fragmentMoonShader]);

      this.vertexTreeShader = initGlShader(gl, "attribute vec2 vPosition;\n\n      varying vec4 fColor;\n\n      void main() {\n        gl_PointSize = 1.0;\n        gl_Position = vec4(vPosition, 0.0, 1.0);\n\n        //fColor = vec4(1.0, 1.0, 1.0, 0.8);\n        fColor = vec4(0.8156862745098039, 0.8470588235294118, 0.8313725490196079, 1.0);\n        //fColor = vec4(0.5490196078431373, 0.5882352941176471, 0.5686274509803921, 0.9);\n      }", gl.VERTEX_SHADER);

      this.fragmentTreeShader = initGlShader(gl, "precision mediump float;\n\n      varying vec4 fColor;\n\n      void main() {\n        gl_FragColor = fColor;\n      }", gl.FRAGMENT_SHADER);

      this.treeProgram = initGlProgram(gl, [this.vertexTreeShader, this.fragmentTreeShader]);
    }
  }, {
    key: "renderGlTree",
    value: function renderGlTree(options) {
      var pointsArray = [];
      var colorsArray = [];
      var pos = options.pos;
      var theta = options.theta;
      var cnt = 0;
      var minH = options.minH * 2 / this.ctx.canvas.height;
      var gl = this.gl;

      function branch(root, h, phi, theta) {
        //theta: 分叉角度
        //phi: root枝杈角度

        var left = {
          x: 0 - h * Math.abs(Math.sin(theta)),
          y: h * Math.abs(Math.cos(theta))
        };

        var right = {
          x: h * Math.abs(Math.sin(theta)),
          y: h * Math.abs(Math.cos(theta))
        };

        //rotate
        var x = left.x * Math.cos(-phi) - left.y * Math.sin(-phi);
        var y = left.x * Math.sin(-phi) + left.y * Math.cos(-phi);

        left.x = x + root.x;
        left.y = y + root.y;

        //rotate
        x = right.x * Math.cos(-phi) - right.y * Math.sin(-phi);
        y = right.x * Math.sin(-phi) + right.y * Math.cos(-phi);

        right.x = x + root.x;
        right.y = y + root.y;

        pointsArray.push(vec2(root.x, root.y));
        pointsArray.push(vec2(left.x, left.y));
        pointsArray.push(vec2(root.x, root.y));
        pointsArray.push(vec2(right.x, right.y));
        cnt++;

        if (h >= minH) {
          branch(left, h * 5 / 6, phi - theta, theta);
          branch(right, h * 5 / 6, phi + theta, theta);
        }
      }

      var h = options.iniH * 2 / this.ctx.canvas.height;
      pointsArray.push(vec2(0.0, -1.0));
      pointsArray.push(vec2(0.0, -1.0 + h));
      branch({ x: 0.0, y: -1.0 + h }, h * 5 / 6, 0, theta * 2 / 3);

      var program = this.treeProgram;
      gl.useProgram(program);
      gl.enable(gl.DEPTH_TEST);

      gl.viewport(pos.x - this.ctx.canvas.width / 2, 0, this.ctx.canvas.width, this.ctx.canvas.height);

      var vBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(pointsArray), gl.STATIC_DRAW);

      var vPosition = gl.getAttribLocation(program, "vPosition");
      gl.vertexAttribPointer(vPosition, 2, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(vPosition);

      gl.drawArrays(gl.LINES, 0, pointsArray.length);
    }
  }, {
    key: "renderGlMoon",
    value: function renderGlMoon(options) {

      function triangle(a, b, c) {
        pointsArray.push(a);
        pointsArray.push(b);
        pointsArray.push(c);

        normalsArray.push(a[0], a[1], a[2], 0.0);
        normalsArray.push(b[0], b[1], b[2], 0.0);
        normalsArray.push(c[0], c[1], c[2], 0.0);
      }

      function divideTriangle(a, b, c, count) {
        if (count > 0) {
          var ab = mix(a, b, 0.5);
          var ac = mix(a, c, 0.5);
          var bc = mix(b, c, 0.5);

          ab = normalize(ab, true);
          ac = normalize(ac, true);
          bc = normalize(bc, true);

          divideTriangle(a, ab, ac, count - 1);
          divideTriangle(ab, b, bc, count - 1);
          divideTriangle(bc, c, ac, count - 1);
          divideTriangle(ab, bc, ac, count - 1);
        } else {
          triangle(a, b, c);
        }
      }

      function tetrahedron(a, b, c, d, n) {
        divideTriangle(a, b, c, n);
        divideTriangle(d, c, b, n);
        divideTriangle(a, d, b, n);
        divideTriangle(a, c, d, n);
      }

      var t = options.t;
      var pos = options.pos;
      var r = options.r;
      var numTimesToSubdivide = options.numTimesToSubdivide || 5;

      if (t >= 1) {
        t = 1.0;
      }
      var phi = 0.1;
      var theta = 0.2 + t * 0.5;

      var gl = this.gl;
      var program = this.moonProgram;

      var pointsArray = [];
      var normalsArray = [];
      var colorArray = [];
      var compIdxArray = [];

      var near = -10;
      var far = 10;
      var radius = 1.5;
      var dr = 5.0 * Math.PI / 180.0;

      var left = -3.0;
      var right = 3.0;
      var ytop = 3.0;
      var bottom = -3.0;

      var va = vec4(0.0, 0.0, -1.0, 1);
      var vb = vec4(0.0, 0.942809, 0.333333, 1);
      var vc = vec4(-0.816497, -0.471405, 0.333333, 1);
      var vd = vec4(0.816497, -0.471405, 0.333333, 1);

      var lightPosition = vec4(1.0, 1.0, 1.0, 0.0);
      var lightAmbient = vec4(1.0, 1.0, 1.0, 1.0);
      var lightDiffuse = vec4(1.0, 1.0, 1.0, 1.0);
      var lightSpecular = vec4(1.0, 1.0, 1.0, 1.0);

      var materialAmbient = vec4(1.0, 1.0, 0.7, 1.0);
      var materialDiffuse = vec4(0.7, 0.7, 0.7, 1.0);
      var materialSpecular = vec4(1.0, 1.0, 1.0, 1.0);
      var materialShininess = 20.0;

      var ambientColor, diffuseColor, specularColor;

      var modelViewMatrix, projectionMatrix;
      var modelViewMatrixLoc, projectionMatrixLoc;

      var normalMatrix, normalMatrixLoc;

      var eye;
      var at = vec3(0.0, 0.0, 0.0);
      var up = vec3(0.0, 1.0, 0.0);

      gl.useProgram(program);
      gl.enable(gl.DEPTH_TEST);

      gl.viewport(pos.x - 2 * r, this.ctx.canvas.height - (pos.y + 2 * r), 4 * r, 4 * r);

      var ambientProduct = mult(lightAmbient, materialAmbient);
      var diffuseProduct = mult(lightDiffuse, materialDiffuse);
      var specularProduct = mult(lightSpecular, materialSpecular);
      var data;

      tetrahedron(va, vb, vc, vd, numTimesToSubdivide);

      var nBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, nBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(normalsArray), gl.STATIC_DRAW);

      var vNormal = gl.getAttribLocation(program, "vNormal");
      gl.vertexAttribPointer(vNormal, 4, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(vNormal);

      var vBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, flatten(pointsArray), gl.STATIC_DRAW);

      var vPosition = gl.getAttribLocation(program, "vPosition");
      gl.vertexAttribPointer(vPosition, 4, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(vPosition);

      modelViewMatrixLoc = gl.getUniformLocation(program, "modelViewMatrix");
      projectionMatrixLoc = gl.getUniformLocation(program, "projectionMatrix");
      normalMatrixLoc = gl.getUniformLocation(program, "normalMatrix");

      gl.uniform4fv(gl.getUniformLocation(program, "ambientProduct"), flatten(ambientProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "diffuseProduct"), flatten(diffuseProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "specularProduct"), flatten(specularProduct));
      gl.uniform4fv(gl.getUniformLocation(program, "lightPosition"), flatten(lightPosition));
      gl.uniform1f(gl.getUniformLocation(program, "shininess"), materialShininess);
      gl.uniform1f(gl.getUniformLocation(program, "t"), t);

      eye = vec3(radius * Math.sin(theta) * Math.cos(phi), radius * Math.sin(theta) * Math.sin(phi), radius * Math.cos(theta));

      modelViewMatrix = lookAt(eye, at, up);
      projectionMatrix = ortho(left, right, bottom, ytop, near, far);

      normalMatrix = [vec3(modelViewMatrix[0][0], modelViewMatrix[0][1], modelViewMatrix[0][2]), vec3(modelViewMatrix[1][0], modelViewMatrix[1][1], modelViewMatrix[1][2]), vec3(modelViewMatrix[2][0], modelViewMatrix[2][1], modelViewMatrix[2][2])];

      gl.uniformMatrix4fv(modelViewMatrixLoc, false, flatten(modelViewMatrix));
      gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix));
      gl.uniformMatrix3fv(normalMatrixLoc, false, flatten(normalMatrix));

      gl.drawArrays(gl.TRIANGLES, 0, pointsArray.length);
    }
  }, {
    key: "clear",
    value: function clear() {
      this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
      this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
    }
  }]);

  return TreeBgSprite;
}(BgSprite);

var Data = function () {
  function Data(editor, text) {
    _classCallCheck(this, Data);

    this._content = [];
    this.editor = editor;
    this._lines = [];
    text && this.addContent(0, text);
  }

  _createClass(Data, [{
    key: "addContent",
    value: function addContent(index, text) {
      window._text = text;
      for (var i = 0; i < text.length; i++) {
        this._content.splice(index, 0, text.charAt(i));
        index++;
      }
      this._lines = this._content.join('').split('\n');
    }
  }, {
    key: "deleteContent",
    value: function deleteContent(index, length) {
      var text = this._content.splice(index - length, length);
      this._lines = this._content.join('').split('\n');
      return text.join('');
    }
  }, {
    key: "getContent",
    value: function getContent(index, length) {
      if (arguments.length === 0) {
        return this._content.join();
      } else {
        return this._content.slice(index, index + length);
      }
    }
  }, {
    key: "getLength",
    value: function getLength() {
      return this._content.length;
    }
  }, {
    key: "getLines",
    value: function getLines() {
      return this._lines;
    }
  }, {
    key: "getLine",
    value: function getLine(index) {
      return this._lines[index];
    }
  }]);

  return Data;
}();

var Editor = function () {
  function Editor(options) {
    _classCallCheck(this, Editor);

    this.options = {
      paddingLeft: 120,
      paddingRight: 120,
      paddingTop: 65,
      paddingBottom: 65,
      bottomBarHeight: 68
    };
    this.onFocus = false;
    this.isToolBarShown = false;
    this.toolBarWidth = 200;
    this.operations = [];
    this.init();
  }

  _createClass(Editor, [{
    key: "init",
    value: function init() {
      this.createStyle();
      this.initData();
      this.selection = new Selection(this);
      this.createElements();
      this.bindEvents();
      this.createMap();
      this.style.initSprites();
    }
  }, {
    key: "resize",
    value: function resize() {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var $editor = $(this.editorNode);

      var width = windowWidth > 800 ? windowWidth - this.options.paddingLeft - this.options.paddingRight : windowWidth - 50;
      var height = windowWidth > 800 ? windowHeight - this.options.bottomBarHeight - this.options.paddingTop - this.options.paddingBottom : windowHeight - 135;
      $editor.css('width', width + 'px');
      $editor.css('height', height + 'px');
      $editor.css('top', (windowHeight - height - this.options.bottomBarHeight) / 2 + 'px');
      this.canvasNode.width = width;
      this.canvasNode.height = height;
      if (this.isToolBarShown && windowWidth > 800) {
        this.updateRelativeLocation((windowHeight - height - this.options.bottomBarHeight) / 2, this.toolBarWidth + (windowWidth - this.toolBarWidth - width) / 2);
      } else {
        this.updateRelativeLocation((windowHeight - height - 75) / 2, (windowWidth - width) / 2);
      }
      this.renderText();
    }
  }, {
    key: "showToolBar",
    value: function showToolBar() {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var width = windowWidth > 800 ? windowWidth - this.options.paddingLeft - this.options.paddingRight : windowWidth - 50;
      var height = windowWidth > 800 ? windowHeight - this.options.bottomBarHeight - this.options.paddingTop - this.options.paddingBottom : windowHeight - 135;
      this.isToolBarShown = true;
      this.updateRelativeLocation((windowHeight - height - this.options.bottomBarHeight) / 2, this.toolBarWidth + (windowWidth - this.toolBarWidth - width) / 2);
    }
  }, {
    key: "hideToolBar",
    value: function hideToolBar() {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var width = windowWidth > 800 ? windowWidth - this.options.paddingLeft - this.options.paddingRight : windowWidth - 50;
      var height = windowWidth > 800 ? windowHeight - this.options.bottomBarHeight - this.options.paddingTop - this.options.paddingBottom : windowHeight - 135;
      this.isToolBarShown = false;
      this.updateRelativeLocation((windowHeight - height - this.options.bottomBarHeight) / 2, (windowWidth - width) / 2);
    }
  }, {
    key: "updateRelativeLocation",
    value: function updateRelativeLocation(top, left) {
      this.top = top;
      this.left = left;
    }
  }, {
    key: "createStyle",
    value: function createStyle() {
      this.style = new Style(this);
    }
  }, {
    key: "createElements",
    value: function createElements() {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var width = windowWidth > 800 ? windowWidth - this.options.paddingLeft - this.options.paddingRight : windowWidth - 50;
      var height = windowWidth > 800 ? windowHeight - this.options.bottomBarHeight - this.options.paddingTop - this.options.paddingBottom : windowHeight - 135;

      this.$editor = $('#editor');
      this.editorNode = this.$editor.get(0);

      this.$editor.css('overflow', 'hidden');
      this.$editor.css('position', 'relative');
      this.$editor.css('box-shadow', '1px 1px 3px 1px #bbb');
      this.$editor.css('margin', 'auto');
      this.$editor.css('width', width + 'px');
      this.$editor.css('height', height + 'px');
      this.$editor.css('top', (windowHeight - height - this.options.bottomBarHeight) / 2 + 'px');
      //this.$editor.css('background-color', 'rgba(217, 218, 214)');
      this.$editor.css('background-color', 'rgb(241, 240, 236)');

      this.hiddenCanvas = document.createElement('canvas');
      this.hiddenCtx = this.hiddenCanvas.getContext('2d');
      this.hiddenCanvas.width = this.style.lineHeight();
      this.hiddenCanvas.height = this.style.lineHeight();

      //for background canvas
      this.$bgCanvas = $('#bg-canvas');
      this.$bgGlcanvas = $('#bg-glcanvas');
      this.bgGlcanvasNode = this.$bgGlcanvas.get(0);
      this.bgCanvasNode = this.$bgCanvas.get(0);
      this.bgCtx = this.bgCanvasNode.getContext('2d');
      this.bgCanvasNode.width = width;
      this.bgCanvasNode.height = height;
      this.bgGlcanvasNode.width = width;
      this.bgGlcanvasNode.height = height;
      this.bgGl = this.bgGlcanvasNode.getContext('webgl', { preserveDrawingBuffer: true });

      //for text canvas
      this.$canvas = $('#text-canvas');
      this.$glcanvas = $('#text-glcanvas');
      this.glcanvasNode = this.$glcanvas.get(0);
      this.canvasNode = this.$canvas.get(0);
      this.ctx = this.canvasNode.getContext('2d');
      this.canvasNode.width = width;
      this.canvasNode.height = height;
      this.glcanvasNode.width = width;
      this.glcanvasNode.height = height;
      this.gl = this.glcanvasNode.getContext('webgl', { preserveDrawingBuffer: true });

      this.$input = $('<textarea>');
      this.inputNode = this.$input.get(0);

      this.$input.css('position', 'absolute');
      this.$input.css('top', '-10px');
      this.$input.css('left', '-10px');
      this.$input.css('width', '1px');
      this.$input.css('height', '1px');

      this.$editor.append(this.$input);

      this.$cursor = $('<div class="cursor"></div>');
      this.cursorNode = this.$cursor.get(0);
      this.$cursor.css('width', '1px');
      this.$cursor.css('height', this.style.lineHeight() + 'px');
      this.$cursor.css('position', 'absolute');
      this.$cursor.css('top', this.selection.rowIndex * this.style.lineHeight());
      this.$cursor.css('left', this.selection.colIndex * this.fontSize);
      this.$cursor.css('background-color', 'black');
      this.$editor.append(this.$cursor);
      this.hideCursor();

      this.updateRelativeLocation((windowHeight - height - this.options.bottomBarHeight) / 2, (windowWidth - width) / 2);
    }
  }, {
    key: "updateCursor",
    value: function updateCursor() {
      var pos = this.selection.getSelEndPosition();
      this.$cursor.css('height', this.style.lineHeight() + 'px');
      this.$cursor.css('left', this.map[pos.rowIndex][pos.colIndex].cursorX + 'px');
      this.$cursor.css('top', this.map[pos.rowIndex][pos.colIndex].cursorY + 'px');
    }
  }, {
    key: "initData",
    value: function initData() {
      this.data = new Data(this);
    }
  }, {
    key: "bindEvents",
    value: function bindEvents() {
      this.$input.on('change', this.onInputChange.bind(this));
      this.$input.on('focus', this.onInputFocus.bind(this));
      this.$input.on('blur', this.onInputBlur.bind(this));
      this.$input.on('keydown', this.onKeyDown.bind(this));
      this.$input.on('compositionstart', this.onCompStart.bind(this));
      this.$input.on('compositionend', this.onCompEnd.bind(this));
      this.$input.on('input', this.onInputChar.bind(this));
      this.$editor.on('touchstart', this.onTouchStart.bind(this));
      this.$editor.on('touchend', this.onTouchEnd.bind(this));
      this.$editor.on('mousedown', this.onMouseDown.bind(this));
      this.$editor.on('mouseup', this.onMouseUp.bind(this));
      $('body').on('click touchend', this.onGlobalClick.bind(this));
      $(window).on('resize', this.resize.bind(this));
      $(window).on('scroll', this.onScroll.bind(this));
    }
  }, {
    key: "findPosfromMap",
    value: function findPosfromMap(x, y) {
      if (this.map.length === 0) {
        return { col: 0, row: 0 };
      }

      var map = this.map;
      var startRow = 0;
      var endRow = this.map.length - 1;
      var row, col;
      var charLength = this.style.fontSize;
      var lineHeight = this.style.lineHeight();
      var pos = { row: 0, col: 0 };

      for (row = 0; row < map.length; row++) {
        for (col = 0; col < map[row].length; col++) {
          if (Math.abs(x - map[row][col].cursorX) < charLength / 2 && y > map[row][col].cursorY && y < map[row][col].cursorY + lineHeight) {
            pos = { row: row, col: col };
            return pos;
          }
        }
      }

      return pos;
    }
  }, {
    key: "clearInputNodeValue",
    value: function clearInputNodeValue() {
      this.inputNode.value = '';
    }
  }, {
    key: "onKeyDown",
    value: function onKeyDown(e) {
      switch (e.keyCode) {
        case 8:
          if (this.map[0].length <= 1 && this.map.length <= 1) {
            return;
          }
          this.deleteText(this.selection.getSelEndPosition().index, this.selection.getSelLength() || 1);
          this.clearInputNodeValue();
          this.selection.update(this.selection.getSelEndPosition().index - 1);
          this.renderText();
          break;
        case 37:
          this.selection.update(this.selection.getSelEndPosition().index - 1);
          this.updateCursor();
          break;
        case 38:

          break;
        case 39:
          this.selection.update(this.selection.getSelEndPosition().index + 1);
          this.updateCursor();
          break;
        case 40:

          break;
      }
    }
  }, {
    key: "insertText",
    value: function insertText(index, text) {
      this.data.addContent(index, text);
      this.selection.update(index + text.length);
      this.operations.push({
        action: index >= this.data.getLength() ? 'add' : 'insert',
        index: index,
        text: text
      });
    }
  }, {
    key: "deleteText",
    value: function deleteText(index, length) {
      var text = this.data.deleteContent(index, length);
      this.operations.push({
        action: 'delete',
        index: index,
        text: text
      });
    }
  }, {
    key: "clearText",
    value: function clearText() {
      this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
      this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
    }
  }, {
    key: "clearCanvas",
    value: function clearCanvas() {
      this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
      this.bgCtx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
      this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
      this.bgGl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
    }
  }, {
    key: "createMap",
    value: function createMap() {
      this.map = [[{
        char: '',
        x: this.style.padding[3],
        y: this.style.padding[0],
        cursorX: this.style.textAlign === 'center' ? this.canvasNode.width / 2 : this.style.padding[3],
        cursorY: this.style.padding[0]
      }]];
    }
  }, {
    key: "_convertWindowPosToCanvas",
    value: function _convertWindowPosToCanvas(x, y) {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var x2 = x - this.left;
      var y2 = y - this.top;

      if (x2 < 0 || x2 > this.canvasNode.width || y2 < 0 || y2 > this.canvasNode.height) {
        x2 = -1;
        y2 = -1;
      }

      return {
        x: x2,
        y: y2
      };
    }
  }, {
    key: "_convertCanvasPosToWindow",
    value: function _convertCanvasPosToWindow(x, y) {
      var windowWidth = $(window).width();
      var windowHeight = $(window).height();
      var x2 = x + this.left;
      var y2 = y + this.top;

      return {
        x: x2,
        y: y2
      };
    }
  }, {
    key: "onGlobalClick",
    value: function onGlobalClick(e) {
      var pos = this._convertWindowPosToCanvas(e.clientX, e.clientY);

      if (pos.x === -1 && pos.y === -1) {
        this.blur();
      }
    }
  }, {
    key: "onInputChange",
    value: function onInputChange(e) {}
  }, {
    key: "onInputBlur",
    value: function onInputBlur(e) {
      this.inputStatus = 'NULL';
    }
  }, {
    key: "onInputFocus",
    value: function onInputFocus(e) {}
  }, {
    key: "onCompStart",
    value: function onCompStart(e) {
      this.inputStatus = 'CHINESE_TYPING';
    }
  }, {
    key: "onCompEnd",
    value: function onCompEnd(e) {
      var that = this;
      setTimeout(function () {
        that.input();
        that.inputStatus = 'CHINESE_TYPE_END';
      }, 100);
    }
  }, {
    key: "onInputChar",
    value: function onInputChar(e) {
      if (this.inputStatus === 'CHINESE_TYPING') {
        return;
      }

      this.inputStatus = 'CHAR_TYPING';
      this.input();
    }
  }, {
    key: "onScroll",
    value: function onScroll(e) {
      this.touch = {};
      this.click = {};
    }
  }, {
    key: "input",
    value: function input(e) {
      this.insertText(this.selection.getSelEndPosition().index, this.inputNode.value);
      this.clearInputNodeValue();
      this.renderText();
    }
  }, {
    key: "onMouseUp",
    value: function onMouseUp(e) {
      e.preventDefault();

      var that = this;

      if (this.click) {
        var pos = this._convertWindowPosToCanvas(e.clientX, e.clientY);
        if (pos.x !== -1 && pos.y !== -1) {
          this.focus(pos.x, pos.y);
        } else {
          this.blur();
        }
      }
      this.click = {};
    }
  }, {
    key: "onMouseDown",
    value: function onMouseDown(e) {
      this.click = {
        x: e.clientX,
        y: e.clientY
      };
    }
  }, {
    key: "onTouchStart",
    value: function onTouchStart(e) {
      this.touch = e.touches[0];
    }
  }, {
    key: "onTouchEnd",
    value: function onTouchEnd(e) {
      e.preventDefault();

      var that = this;

      if (this.touch) {
        var pos = this._convertWindowPosToCanvas(e.changedTouches[0].clientX, e.changedTouches[0].clientY);
        if (pos.x !== -1 && pos.y !== -1) {
          this.focus(pos.x, pos.y);
        } else {
          this.blur();
        }
      }
      this.touch = {};
    }
  }, {
    key: "hideCursor",
    value: function hideCursor() {
      this.$cursor.css('visibility', 'hidden');
    }
  }, {
    key: "blur",
    value: function blur() {
      this.hideCursor();
      this.inputNode.blur();
      this.onFocus = false;
    }
  }, {
    key: "focus",
    value: function focus(x, y) {
      var pos = this.findPosfromMap(x, y);
      this.selection.update(pos.row, pos.col);
      this.updateCursor();
      this.$input.focus();
      this.$cursor.css('visibility', 'visible');
      this.onFocus = true;
    }
  }, {
    key: "updateMap",
    value: function updateMap(map) {
      this.map = map;
    }
  }, {
    key: "measureStrLength",
    value: function measureStrLength(str) {
      var length = 0;
      var chars = str.split('');
      chars.forEach(function (char) {
        if (char) {
          var code = char.charCodeAt(0);
          if (code >= 48 && code <= 57 || code >= 65 && code <= 90 || code >= 97 && code <= 122) {
            length += 0.5;
          } else {
            length += 1;
          }
        }
      });

      return length;
    }
  }, {
    key: "renderText",
    value: function renderText() {
      this.bgSprite.drawStatic();

      var xCursor = this.style.padding[3];
      var yCursor = this.style.padding[0];
      var map = [];
      var lines = this.data.getLines();

      if (lines.length === 0) {
        lines.push('');
      }

      this.ctx.font = this.style.fontSize + "px serif";
      this.ctx.textAlign = 'left';
      this.ctx.textBaseline = 'middle';

      switch (this.style.textAlign) {
        case 'left':
          for (var i = 0; i < lines.length; i++) {
            var lineMap = [];
            var charLength = this.style.fontSize;
            xCursor = this.style.padding[3];

            var chars = lines[i].split('');
            lineMap.push({
              char: '',
              x: xCursor,
              y: yCursor,
              cursorX: xCursor,
              cursorY: yCursor
            });
            for (var j = 0; j < chars.length; j++) {
              var charGap = charLength;
              xCursor += charLength;
              lineMap.push({
                char: chars[j],
                x: xCursor - charGap,
                y: yCursor,
                cursorX: xCursor,
                cursorY: yCursor
              });
              if (j + 1 < chars.length && xCursor + charGap > this.canvasNode.width - this.style.padding[1] - this.style.padding[3]) {
                yCursor += this.style.lineHeight();
                xCursor = this.style.padding[3];
              }
            }
            map.push(lineMap);
            yCursor += this.style.lineHeight();
            yCursor += this.style.space;
          }
          break;
        case 'center':
          for (var i = 0; i < lines.length; i++) {
            var lineMap = [];
            var chars = lines[i].split('');
            var charLength = this.style.fontSize; //this.ctx.measureText(chars[j]).width;
            var length = charLength * chars.length; //this.measureStrLength(lines[i]);//this.ctx.measureText(lines[i]).width;

            if (length > this.canvasNode.width - this.style.padding[0] - this.style.padding[3]) {
              xCursor = this.style.padding[3];
            } else {
              xCursor = (this.canvasNode.width - this.style.padding[0] - this.style.padding[3] - length) / 2 + this.style.padding[3];
            }
            lineMap.push({
              char: '',
              x: xCursor,
              y: yCursor,
              cursorX: xCursor,
              cursorY: yCursor
            });

            while (length > this.canvasNode.width - this.style.padding[0] - this.style.padding[3]) {
              xCursor = this.style.padding[3];
              for (var j = 0; j < chars.length; j++) {
                var charGap = charLength; // * this.measureStrLength(chars[j]);
                xCursor += charGap;
                lineMap.push({
                  char: chars[j],
                  x: xCursor - charGap,
                  y: yCursor,
                  cursorX: xCursor,
                  cursorY: yCursor
                });
                if (j + 1 < chars.length && xCursor + charGap > this.canvasNode.width - this.style.padding[1]) {
                  yCursor += this.style.lineHeight();
                  chars.splice(0, j + 1);
                  length = charLength * chars.length; //maersureStrLength(chars.join(''));
                  break;
                }
              }
            }

            if (length) {
              xCursor = (this.canvasNode.width - this.style.padding[0] - this.style.padding[3] - length) / 2 + this.style.padding[3];
              for (j = 0; j < chars.length; j++) {
                var charGap = charLength; // * this.measureStrLength(chars[j]);
                xCursor += charGap;
                lineMap.push({
                  char: chars[j],
                  x: xCursor - charGap,
                  y: yCursor,
                  cursorX: xCursor,
                  cursorY: yCursor
                });
              }
            }

            map.push(lineMap);
            yCursor += this.style.lineHeight();
            yCursor += this.style.space;
          }
          break;
      }
      this.updateMap(map);
      this.textSprite.update();
      this._fillText();
      this.updateCursor();
    }
  }, {
    key: "_fillText",
    value: function _fillText() {
      if (this.map.length === 1 && this.map[0].length === 1) {
        this.textSprite.clear();
      } else {
        $('.render-tip').addClass('show');
        setTimeout(this.textSprite.drawStatic.bind(this.textSprite), 10);
      }
    }
  }, {
    key: "play",
    value: function play() {
      this.animating = true;
      this.animationInfo = {
        textStop: false,
        bgStop: false
      };
      this.startTime = Date.now();
      this.textSprite.update();
      this.bgSprite.update();

      window.requestAnimationFrame(this.tick.bind(this));
    }
  }, {
    key: "tick",
    value: function tick() {
      if (!this.animating) {
        return;
      }

      var t = Date.now() - this.startTime;
      !this.animationInfo.textStop && (this.animationInfo.textStop = this.textSprite.advance(t));
      !this.animationInfo.bgStop && (this.animationInfo.bgStop = this.bgSprite.advance(t));

      if (this.animationInfo.textStop && this.animationInfo.bgStop) {
        this.stopPlay();
      } else {
        this.animationInfo.bgStop ? this.bgSprite.drawStatic() : this.bgSprite.drawFrame();
        this.animationInfo.textStop ? this.textSprite.drawStatic() : this.textSprite.drawFrame();
        window.requestAnimationFrame(this.tick.bind(this));
      }
    }
  }, {
    key: "stopPlay",
    value: function stopPlay() {
      this.animating = false;
      this.renderText();
      bottomBar.config.$dom.find('.play').removeClass('stop');
    }
  }, {
    key: "generateGif",
    value: function generateGif() {}
  }, {
    key: "generateJpeg",
    value: function generateJpeg() {
      var canvas = document.createElement('canvas');
      canvas.width = this.canvasNode.width;
      canvas.height = this.canvasNode.height;
      var ctx = canvas.getContext('2d');
      ctx.drawImage(this.bgCanvasNode, 0, 0);
      ctx.drawImage(this.bgGlcanvasNode, 0, 0);
      ctx.drawImage(this.canvasNode, 0, 0);
      ctx.drawImage(this.glcanvasNode, 0, 0);

      var imgData = canvas.toDataURL("image/jpeg", 1.0);
      return imgData;
    }
  }, {
    key: "generatePng",
    value: function generatePng() {
      var canvas = document.createElement('canvas');
      canvas.width = this.canvasNode.width;
      canvas.height = this.canvasNode.height;
      var ctx = canvas.getContext('2d');
      ctx.drawImage(this.bgCanvasNode, 0, 0);
      ctx.drawImage(this.bgGlcanvasNode, 0, 0);
      ctx.drawImage(this.canvasNode, 0, 0);
      ctx.drawImage(this.glcanvasNode, 0, 0);

      var imgData = canvas.toDataURL("image/png");
      return imgData;
    }
  }]);

  return Editor;
}();

var Selection = function () {
  function Selection(editor) {
    _classCallCheck(this, Selection);

    this.editor = editor;
    this.startColIndex = this.endColIndex = 0;
    this.startRowIndex = this.endRowIndex = 0;
  }

  _createClass(Selection, [{
    key: "getSelEndPosition",
    value: function getSelEndPosition() {
      return {
        colIndex: this.endColIndex,
        rowIndex: this.endRowIndex,
        index: this._convertIndexByMap(this.endRowIndex, this.endColIndex)
      };
    }
  }, {
    key: "getSelLength",
    value: function getSelLength() {
      return 0;
    }
  }, {
    key: "update",
    value: function update(startRowIndex, startColIndex, endRowIndex, endColIndex) {
      if (arguments.length === 1) {
        var index = arguments[0];
        if (index > this.editor.data.getLength() || index < 0) {
          return;
        }
        this.index = index;
        var pos = this._convertMapByIndex(index);
        this.startColIndex = this.endColIndex = pos.col;
        this.startRowIndex = this.endRowIndex = pos.row;
      } else if (arguments.length === 2) {
        this.startColIndex = this.endColIndex = startColIndex;
        this.startRowIndex = this.endRowIndex = startRowIndex;
      }
    }
  }, {
    key: "_convertIndexByMap",
    value: function _convertIndexByMap(row, col) {
      var lines = this.editor.data.getLines();
      var index = 0;
      for (var i = 0; i < lines.length; i++) {
        if (i < row) {
          index += lines[i].length + 1;
        } else if (i === row) {
          index += col;
        }
      }
      return index;
    }
  }, {
    key: "_convertMapByIndex",
    value: function _convertMapByIndex(index) {
      var lines = this.editor.data.getLines();
      var row = 0;
      var col = 0;

      for (var i = 0; i < lines.length; i++) {
        if (index > lines[i].length) {
          index -= lines[i].length + 1;
          row++;
        } else {
          col = index;
          break;
        }
      }
      return {
        row: row,
        col: col
      };
    }
  }]);

  return Selection;
}();

var Style = function () {
  function Style(editor) {
    _classCallCheck(this, Style);

    this.editor = editor;
    this.initDefaultStyle();
  }

  _createClass(Style, [{
    key: "changeStyle",
    value: function changeStyle(options) {
      var keys = Object.keys(options);
      for (var i = 0; i < keys.length; i++) {
        this[keys[i]] = options[keys[i]];

        if (keys[i] === 'textStyle') {
          this.changeTextStyle(options[keys[i]]);
        }

        if (keys[i] === 'background') {
          this.changeBgStyle(options[keys[i]]);
        }
      }

      if (this.editor.hiddenCanvas) {
        this.editor.hiddenCanvas.width = this.lineHeight();
        this.editor.hiddenCanvas.height = this.lineHeight();
      }
    }
  }, {
    key: "changeTextStyle",
    value: function changeTextStyle(index) {
      this.editor.textSprite = util.getSpriteEntity(config.textStyleMap[index].Klass, this.editor);
      this.editor.textSprite.changeTextStyle(config.textStyleMap[index].style);
      config.animationMap = this.editor.textSprite.getAnimationMap();
      this.editor.textSprite.update();
    }
  }, {
    key: "changeBgStyle",
    value: function changeBgStyle(index) {
      this.editor.bgSprite = util.getSpriteEntity(config.backgroundMap[index].Klass, this.editor);
      this.editor.bgSprite.update();
    }
  }, {
    key: "lineHeight",
    value: function lineHeight() {
      return this.lineHeightRatio * this.fontSize;
    }
  }, {
    key: "initDefaultStyle",
    value: function initDefaultStyle() {
      this.changeStyle({
        textAlign: util.getInitialState('textAlign'),
        fontSize: util.getInitialState('fontSize'),
        lineHeightRatio: 1.2,
        textWeight: 800,
        padding: [80, 50, 50, 80],
        space: 20,
        animation: util.getInitialState('animation'),
        fontColor: util.getInitialState('fontColor')
      });
    }
  }, {
    key: "initSprites",
    value: function initSprites() {
      this.changeStyle({
        textStyle: util.getInitialState('textStyle'),
        background: util.getInitialState('background')
      });
    }
  }]);

  return Style;
}();

var bottomBar = {
  config: {
    $dom: $('.bottom-bar')
  },

  init: function init() {
    this.bindEvents();
  },

  bindEvents: function bindEvents() {
    var that = this;

    that.config.$dom.find('.export').on('click touchend', function (e) {
      selectModal.show();
    });

    that.config.$dom.find('.play').on('touchend click', function (e) {
      if ($(this).hasClass('stop')) {
        editor.stopPlay();
        $(this).removeClass('stop');
      } else {
        editor.play();
        $(this).addClass('stop');
      }
    });

    that.config.$dom.find('.tool').on('touchend click', function (e) {
      if (toolBar.config.isShown) {
        toolBar.hide(editor);
        $(this).find('.label').text('工具栏');
      } else {
        toolBar.show(editor);
        $(this).find('.label').text('隐藏');
      }
    });
  }
};
var resultModal = {
  config: {
    $dom: $('.export-result-modal')
  },

  init: function init() {
    this.bindEvents();
  },

  bindEvents: function bindEvents() {
    var that = this;
    that.config.$dom.find('.dim-layer').on('click touchend', function (e) {
      that.hide();
    });
  },

  show: function show(imgData) {
    this.config.$dom.find('img').attr('src', imgData);
    this.config.$dom.addClass('show');
  },

  hide: function hide() {
    this.config.$dom.removeClass('show');
  }
};
var selectModal = {
  config: {
    $dom: $('.export-select-modal'),
    selection: 'jpeg'
  },

  init: function init() {
    this.bindEvents();
  },

  bindEvents: function bindEvents() {
    var that = this;

    that.config.$dom.find('.item').on('click touchend', function (e) {
      that.config.selection = $(this).text().toLowerCase();
      that.config.$dom.find('.item').removeClass('active');
      $(this).addClass('active');
    });

    that.config.$dom.find('.confirm').on('click touchend', function (e) {
      that.hide();
      that.export(that.config.selection);
    });

    that.config.$dom.find('.cancel').on('click touchend', function (e) {
      that.hide();
    });
  },

  show: function show() {
    this.config.$dom.addClass('show');
  },

  hide: function hide() {
    this.config.$dom.removeClass('show');
  },

  export: function _export(option) {
    var imgData;
    switch (option) {
      case 'jpeg':
        imgData = editor.generateJpeg();
        resultModal.show(imgData);
        break;
      case 'png':
        imgData = editor.generatePng();
        resultModal.show(imgData);
        break;
      case 'gif':
        editor.generateGif().then(function (imgData) {
          resultModal.show(imgData);
        });
        break;
    }
  }
};
var toolBar = {
  config: {
    isShown: false
  },

  init: function init() {
    this.onChange = false;
    this.config.$dom.find('.tool-bar-item').each(function () {
      var classes = $(this).classes();
      var re = /i-(.*)/i;
      $(this).find('.value').append('<span class="spinner-icon"></span>');
      $(this).find('.spinner-icon').hide();

      classes.forEach(function (className) {
        if (re.test(className)) {
          var name = RegExp.$1;
          toolBar.update(name);
        }
      });
    });
    this.updateBgColor();
    this.bindEvents();
  },

  bindEvents: function bindEvents() {
    var that = this;

    that.config.$dom.find('.tool-bar-item .value').on('click touchend', function (e) {
      var classes = $(this).parent().classes();
      var $el = $(this);
      var re = /i-(.*)/i;

      classes.forEach(function (className) {
        if (re.test(className)) {
          var name = RegExp.$1;
          if (!that[name + 'OnChange']) {
            $el.find('.spinner-icon').show();
            toolBar.changeStyle(name);
          }
        }
      });
    });

    that.config.$dom.find('.close-icon').on('touchend click', function (e) {
      that.hide(editor);
      bottomBar.config.$dom.find('.label').text('工具栏');
    });

    that.config.$dom.find('.color-switch').on('touchend click', function (e) {
      that.changeBgColor();
    });
  },

  hide: function hide(editor) {
    $('.wrapper').removeClass('show-tool-bar');
    $('.wrapper').addClass('hide-tool-bar');
    editor.hideToolBar();
    this.config.isShown = false;
  },

  show: function show(editor) {
    if ($('.wrapper').hasClass('default')) {
      $('.wrapper').removeClass('default');
    } else {
      $('.wrapper').removeClass('hide-tool-bar');
    }
    $('.wrapper').addClass('show-tool-bar');
    editor.showToolBar();
    this.config.isShown = true;
  },

  nextIndex: function nextIndex(name) {
    var index = config.state[name + 'Index'];
    var map = config[name + 'Map'];
    if (index < map.length - 1) {
      config.state[name + 'Index'] += 1;
    } else {
      config.state[name + 'Index'] = 0;
    }
    return config.state[name + 'Index'];
  },

  changeStyle: function changeStyle(name) {
    var index = this.nextIndex(name);
    this.update(name);
  },

  update: function update(name) {
    this[name + 'OnChange'] = true;
    var index = config.state[name + 'Index'];
    var style = config[name + 'Map'][index];
    var value = style.value;
    var label = style.label;
    var options = {};
    options[name] = value;
    $('.i-' + name + ' .value .text').text(label);

    var that = this;
    editor.style.changeStyle(options);

    if (this[name + 'Change'] && typeof this[name + 'Change'] === 'function') {
      this[name + 'Change'].call(this, style);
    }

    editor.clearCanvas();
    editor.renderText();

    setTimeout(function (name, label) {
      that[name + 'OnChange'] = false;
      $('.i-' + name + ' .value .text').text(label);
      $('.i-' + name + ' .value .spinner-icon').hide();
    }, 500, name, label);
  },

  textStyleChange: function textStyleChange(style) {
    config.state['animationIndex'] = 0;
    this.update('animation');
  },

  backgroundChange: function backgroundChange(style) {
    var colors = style.colors;
    if (colors.length === 0) {
      $('.color-switch').hide();
      config.state.bgColorIndex = null;
    } else {
      $('.color-switch').css('background-color', colors[0]).show();
      config.state.bgColorIndex = 0;
    }
  },

  changeBgColor: function changeBgColor() {
    var colors = config.backgroundMap[config.state.backgroundIndex].colors;
    var index = config.state.bgColorIndex;
    index = index + 1 >= colors.length ? 0 : ++index;
    config.state.bgColorIndex = index;
    this.updateBgColor();
  },

  updateBgColor: function updateBgColor() {
    editor.renderText();
    var colors = config.backgroundMap[config.state.backgroundIndex].colors;
    var index = config.state.bgColorIndex;
    $('.color-switch').css('background-color', colors[index]);
  }
};

toolBar.config.$dom = $('.tool-bar-wrap');

(function e(t, n, r) {
  function s(o, u) {
    if (!n[o]) {
      if (!t[o]) {
        var a = typeof require == "function" && require;if (!u && a) return a(o, !0);if (i) return i(o, !0);var f = new Error("Cannot find module '" + o + "'");throw f.code = "MODULE_NOT_FOUND", f;
      }var l = n[o] = { exports: {} };t[o][0].call(l.exports, function (e) {
        var n = t[o][1][e];return s(n ? n : e);
      }, l, l.exports, e, t, n, r);
    }return n[o].exports;
  }var i = typeof require == "function" && require;for (var o = 0; o < r.length; o++) {
    s(r[o]);
  }return s;
})({ 1: [function (require, module, exports) {
    window.earcut = require('earcut');
    window.editor = new Editor({});

    toolBar.init();
    bottomBar.init();
    resultModal.init();
    selectModal.init();
  }, { "earcut": 2 }], 2: [function (require, module, exports) {
    'use strict';

    module.exports = earcut;

    function earcut(data, holeIndices, dim) {

      dim = dim || 2;

      var hasHoles = holeIndices && holeIndices.length,
          outerLen = hasHoles ? holeIndices[0] * dim : data.length,
          outerNode = linkedList(data, 0, outerLen, dim, true),
          triangles = [];

      if (!outerNode) return triangles;

      var minX, minY, maxX, maxY, x, y, size;

      if (hasHoles) outerNode = eliminateHoles(data, holeIndices, outerNode, dim);

      // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
      if (data.length > 80 * dim) {
        minX = maxX = data[0];
        minY = maxY = data[1];

        for (var i = dim; i < outerLen; i += dim) {
          x = data[i];
          y = data[i + 1];
          if (x < minX) minX = x;
          if (y < minY) minY = y;
          if (x > maxX) maxX = x;
          if (y > maxY) maxY = y;
        }

        // minX, minY and size are later used to transform coords into integers for z-order calculation
        size = Math.max(maxX - minX, maxY - minY);
      }

      earcutLinked(outerNode, triangles, dim, minX, minY, size);

      return triangles;
    }

    // create a circular doubly linked list from polygon points in the specified winding order
    function linkedList(data, start, end, dim, clockwise) {
      var i, last;

      if (clockwise === signedArea(data, start, end, dim) > 0) {
        for (i = start; i < end; i += dim) {
          last = insertNode(i, data[i], data[i + 1], last);
        }
      } else {
        for (i = end - dim; i >= start; i -= dim) {
          last = insertNode(i, data[i], data[i + 1], last);
        }
      }

      if (last && equals(last, last.next)) {
        removeNode(last);
        last = last.next;
      }

      return last;
    }

    // eliminate colinear or duplicate points
    function filterPoints(start, end) {
      if (!start) return start;
      if (!end) end = start;

      var p = start,
          again;
      do {
        again = false;

        if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) === 0)) {
          removeNode(p);
          p = end = p.prev;
          if (p === p.next) return null;
          again = true;
        } else {
          p = p.next;
        }
      } while (again || p !== end);

      return end;
    }

    // main ear slicing loop which triangulates a polygon (given as a linked list)
    function earcutLinked(ear, triangles, dim, minX, minY, size, pass) {
      if (!ear) return;

      // interlink polygon nodes in z-order
      if (!pass && size) indexCurve(ear, minX, minY, size);

      var stop = ear,
          prev,
          next;

      // iterate through ears, slicing them one by one
      while (ear.prev !== ear.next) {
        prev = ear.prev;
        next = ear.next;

        if (size ? isEarHashed(ear, minX, minY, size) : isEar(ear)) {
          // cut off the triangle
          triangles.push(prev.i / dim);
          triangles.push(ear.i / dim);
          triangles.push(next.i / dim);

          removeNode(ear);

          // skipping the next vertice leads to less sliver triangles
          ear = next.next;
          stop = next.next;

          continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if (ear === stop) {
          // try filtering points and slicing again
          if (!pass) {
            earcutLinked(filterPoints(ear), triangles, dim, minX, minY, size, 1);

            // if this didn't work, try curing all small self-intersections locally
          } else if (pass === 1) {
            ear = cureLocalIntersections(ear, triangles, dim);
            earcutLinked(ear, triangles, dim, minX, minY, size, 2);

            // as a last resort, try splitting the remaining polygon into two
          } else if (pass === 2) {
            splitEarcut(ear, triangles, dim, minX, minY, size);
          }

          break;
        }
      }
    }

    // check whether a polygon node forms a valid ear with adjacent nodes
    function isEar(ear) {
      var a = ear.prev,
          b = ear,
          c = ear.next;

      if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

      // now make sure we don't have other points inside the potential ear
      var p = ear.next.next;

      while (p !== ear.prev) {
        if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
        p = p.next;
      }

      return true;
    }

    function isEarHashed(ear, minX, minY, size) {
      var a = ear.prev,
          b = ear,
          c = ear.next;

      if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

      // triangle bbox; min & max are calculated like this for speed
      var minTX = a.x < b.x ? a.x < c.x ? a.x : c.x : b.x < c.x ? b.x : c.x,
          minTY = a.y < b.y ? a.y < c.y ? a.y : c.y : b.y < c.y ? b.y : c.y,
          maxTX = a.x > b.x ? a.x > c.x ? a.x : c.x : b.x > c.x ? b.x : c.x,
          maxTY = a.y > b.y ? a.y > c.y ? a.y : c.y : b.y > c.y ? b.y : c.y;

      // z-order range for the current triangle bbox;
      var minZ = zOrder(minTX, minTY, minX, minY, size),
          maxZ = zOrder(maxTX, maxTY, minX, minY, size);

      // first look for points inside the triangle in increasing z-order
      var p = ear.nextZ;

      while (p && p.z <= maxZ) {
        if (p !== ear.prev && p !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
        p = p.nextZ;
      }

      // then look for points in decreasing z-order
      p = ear.prevZ;

      while (p && p.z >= minZ) {
        if (p !== ear.prev && p !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
        p = p.prevZ;
      }

      return true;
    }

    // go through all polygon nodes and cure small local self-intersections
    function cureLocalIntersections(start, triangles, dim) {
      var p = start;
      do {
        var a = p.prev,
            b = p.next.next;

        if (!equals(a, b) && intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {

          triangles.push(a.i / dim);
          triangles.push(p.i / dim);
          triangles.push(b.i / dim);

          // remove two nodes involved
          removeNode(p);
          removeNode(p.next);

          p = start = b;
        }
        p = p.next;
      } while (p !== start);

      return p;
    }

    // try splitting polygon into two and triangulate them independently
    function splitEarcut(start, triangles, dim, minX, minY, size) {
      // look for a valid diagonal that divides the polygon into two
      var a = start;
      do {
        var b = a.next.next;
        while (b !== a.prev) {
          if (a.i !== b.i && isValidDiagonal(a, b)) {
            // split the polygon in two by the diagonal
            var c = splitPolygon(a, b);

            // filter colinear points around the cuts
            a = filterPoints(a, a.next);
            c = filterPoints(c, c.next);

            // run earcut on each half
            earcutLinked(a, triangles, dim, minX, minY, size);
            earcutLinked(c, triangles, dim, minX, minY, size);
            return;
          }
          b = b.next;
        }
        a = a.next;
      } while (a !== start);
    }

    // link every hole into the outer loop, producing a single-ring polygon without holes
    function eliminateHoles(data, holeIndices, outerNode, dim) {
      var queue = [],
          i,
          len,
          start,
          end,
          list;

      for (i = 0, len = holeIndices.length; i < len; i++) {
        start = holeIndices[i] * dim;
        end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
        list = linkedList(data, start, end, dim, false);
        if (list === list.next) list.steiner = true;
        queue.push(getLeftmost(list));
      }

      queue.sort(compareX);

      // process holes from left to right
      for (i = 0; i < queue.length; i++) {
        eliminateHole(queue[i], outerNode);
        outerNode = filterPoints(outerNode, outerNode.next);
      }

      return outerNode;
    }

    function compareX(a, b) {
      return a.x - b.x;
    }

    // find a bridge between vertices that connects hole with an outer ring and and link it
    function eliminateHole(hole, outerNode) {
      outerNode = findHoleBridge(hole, outerNode);
      if (outerNode) {
        var b = splitPolygon(outerNode, hole);
        filterPoints(b, b.next);
      }
    }

    // David Eberly's algorithm for finding a bridge between hole and outer polygon
    function findHoleBridge(hole, outerNode) {
      var p = outerNode,
          hx = hole.x,
          hy = hole.y,
          qx = -Infinity,
          m;

      // find a segment intersected by a ray from the hole's leftmost point to the left;
      // segment's endpoint with lesser x will be potential connection point
      do {
        if (hy <= p.y && hy >= p.next.y) {
          var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
          if (x <= hx && x > qx) {
            qx = x;
            if (x === hx) {
              if (hy === p.y) return p;
              if (hy === p.next.y) return p.next;
            }
            m = p.x < p.next.x ? p : p.next;
          }
        }
        p = p.next;
      } while (p !== outerNode);

      if (!m) return null;

      if (hx === qx) return m.prev; // hole touches outer segment; pick lower endpoint

      // look for points inside the triangle of hole point, segment intersection and endpoint;
      // if there are no points found, we have a valid connection;
      // otherwise choose the point of the minimum angle with the ray as connection point

      var stop = m,
          mx = m.x,
          my = m.y,
          tanMin = Infinity,
          tan;

      p = m.next;

      while (p !== stop) {
        if (hx >= p.x && p.x >= mx && pointInTriangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {

          tan = Math.abs(hy - p.y) / (hx - p.x); // tangential

          if ((tan < tanMin || tan === tanMin && p.x > m.x) && locallyInside(p, hole)) {
            m = p;
            tanMin = tan;
          }
        }

        p = p.next;
      }

      return m;
    }

    // interlink polygon nodes in z-order
    function indexCurve(start, minX, minY, size) {
      var p = start;
      do {
        if (p.z === null) p.z = zOrder(p.x, p.y, minX, minY, size);
        p.prevZ = p.prev;
        p.nextZ = p.next;
        p = p.next;
      } while (p !== start);

      p.prevZ.nextZ = null;
      p.prevZ = null;

      sortLinked(p);
    }

    // Simon Tatham's linked list merge sort algorithm
    // http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
    function sortLinked(list) {
      var i,
          p,
          q,
          e,
          tail,
          numMerges,
          pSize,
          qSize,
          inSize = 1;

      do {
        p = list;
        list = null;
        tail = null;
        numMerges = 0;

        while (p) {
          numMerges++;
          q = p;
          pSize = 0;
          for (i = 0; i < inSize; i++) {
            pSize++;
            q = q.nextZ;
            if (!q) break;
          }

          qSize = inSize;

          while (pSize > 0 || qSize > 0 && q) {

            if (pSize === 0) {
              e = q;
              q = q.nextZ;
              qSize--;
            } else if (qSize === 0 || !q) {
              e = p;
              p = p.nextZ;
              pSize--;
            } else if (p.z <= q.z) {
              e = p;
              p = p.nextZ;
              pSize--;
            } else {
              e = q;
              q = q.nextZ;
              qSize--;
            }

            if (tail) tail.nextZ = e;else list = e;

            e.prevZ = tail;
            tail = e;
          }

          p = q;
        }

        tail.nextZ = null;
        inSize *= 2;
      } while (numMerges > 1);

      return list;
    }

    // z-order of a point given coords and size of the data bounding box
    function zOrder(x, y, minX, minY, size) {
      // coords are transformed into non-negative 15-bit integer range
      x = 32767 * (x - minX) / size;
      y = 32767 * (y - minY) / size;

      x = (x | x << 8) & 0x00FF00FF;
      x = (x | x << 4) & 0x0F0F0F0F;
      x = (x | x << 2) & 0x33333333;
      x = (x | x << 1) & 0x55555555;

      y = (y | y << 8) & 0x00FF00FF;
      y = (y | y << 4) & 0x0F0F0F0F;
      y = (y | y << 2) & 0x33333333;
      y = (y | y << 1) & 0x55555555;

      return x | y << 1;
    }

    // find the leftmost node of a polygon ring
    function getLeftmost(start) {
      var p = start,
          leftmost = start;
      do {
        if (p.x < leftmost.x) leftmost = p;
        p = p.next;
      } while (p !== start);

      return leftmost;
    }

    // check if a point lies within a convex triangle
    function pointInTriangle(ax, ay, bx, by, cx, cy, px, py) {
      return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 && (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 && (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
    }

    // check if a diagonal between two polygon nodes is valid (lies in polygon interior)
    function isValidDiagonal(a, b) {
      return a.next.i !== b.i && a.prev.i !== b.i && !intersectsPolygon(a, b) && locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b);
    }

    // signed area of a triangle
    function area(p, q, r) {
      return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    }

    // check if two points are equal
    function equals(p1, p2) {
      return p1.x === p2.x && p1.y === p2.y;
    }

    // check if two segments intersect
    function intersects(p1, q1, p2, q2) {
      if (equals(p1, q1) && equals(p2, q2) || equals(p1, q2) && equals(p2, q1)) return true;
      return area(p1, q1, p2) > 0 !== area(p1, q1, q2) > 0 && area(p2, q2, p1) > 0 !== area(p2, q2, q1) > 0;
    }

    // check if a polygon diagonal intersects any polygon segments
    function intersectsPolygon(a, b) {
      var p = a;
      do {
        if (p.i !== a.i && p.next.i !== a.i && p.i !== b.i && p.next.i !== b.i && intersects(p, p.next, a, b)) return true;
        p = p.next;
      } while (p !== a);

      return false;
    }

    // check if a polygon diagonal is locally inside the polygon
    function locallyInside(a, b) {
      return area(a.prev, a, a.next) < 0 ? area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 : area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
    }

    // check if the middle point of a polygon diagonal is inside the polygon
    function middleInside(a, b) {
      var p = a,
          inside = false,
          px = (a.x + b.x) / 2,
          py = (a.y + b.y) / 2;
      do {
        if (p.y > py !== p.next.y > py && px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x) inside = !inside;
        p = p.next;
      } while (p !== a);

      return inside;
    }

    // link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
    // if one belongs to the outer ring and another to a hole, it merges it into a single ring
    function splitPolygon(a, b) {
      var a2 = new Node(a.i, a.x, a.y),
          b2 = new Node(b.i, b.x, b.y),
          an = a.next,
          bp = b.prev;

      a.next = b;
      b.prev = a;

      a2.next = an;
      an.prev = a2;

      b2.next = a2;
      a2.prev = b2;

      bp.next = b2;
      b2.prev = bp;

      return b2;
    }

    // create a node and optionally link it with previous one (in a circular doubly linked list)
    function insertNode(i, x, y, last) {
      var p = new Node(i, x, y);

      if (!last) {
        p.prev = p;
        p.next = p;
      } else {
        p.next = last.next;
        p.prev = last;
        last.next.prev = p;
        last.next = p;
      }
      return p;
    }

    function removeNode(p) {
      p.next.prev = p.prev;
      p.prev.next = p.next;

      if (p.prevZ) p.prevZ.nextZ = p.nextZ;
      if (p.nextZ) p.nextZ.prevZ = p.prevZ;
    }

    function Node(i, x, y) {
      // vertice index in coordinates array
      this.i = i;

      // vertex coordinates
      this.x = x;
      this.y = y;

      // previous and next vertice nodes in a polygon ring
      this.prev = null;
      this.next = null;

      // z-order curve value
      this.z = null;

      // previous and next nodes in z-order
      this.prevZ = null;
      this.nextZ = null;

      // indicates whether this is a steiner point
      this.steiner = false;
    }

    // return a percentage difference between the polygon area and its triangulation area;
    // used to verify correctness of triangulation
    earcut.deviation = function (data, holeIndices, dim, triangles) {
      var hasHoles = holeIndices && holeIndices.length;
      var outerLen = hasHoles ? holeIndices[0] * dim : data.length;

      var polygonArea = Math.abs(signedArea(data, 0, outerLen, dim));
      if (hasHoles) {
        for (var i = 0, len = holeIndices.length; i < len; i++) {
          var start = holeIndices[i] * dim;
          var end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
          polygonArea -= Math.abs(signedArea(data, start, end, dim));
        }
      }

      var trianglesArea = 0;
      for (i = 0; i < triangles.length; i += 3) {
        var a = triangles[i] * dim;
        var b = triangles[i + 1] * dim;
        var c = triangles[i + 2] * dim;
        trianglesArea += Math.abs((data[a] - data[c]) * (data[b + 1] - data[a + 1]) - (data[a] - data[b]) * (data[c + 1] - data[a + 1]));
      }

      return polygonArea === 0 && trianglesArea === 0 ? 0 : Math.abs((trianglesArea - polygonArea) / polygonArea);
    };

    function signedArea(data, start, end, dim) {
      var sum = 0;
      for (var i = start, j = end - dim; i < end; i += dim) {
        sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
        j = i;
      }
      return sum;
    }

    // turn a polygon in a multi-dimensional array form (e.g. as in GeoJSON) into a form Earcut accepts
    earcut.flatten = function (data) {
      var dim = data[0][0].length,
          result = { vertices: [], holes: [], dimensions: dim },
          holeIndex = 0;

      for (var i = 0; i < data.length; i++) {
        for (var j = 0; j < data[i].length; j++) {
          for (var d = 0; d < dim; d++) {
            result.vertices.push(data[i][j][d]);
          }
        }
        if (i > 0) {
          holeIndex += data[i - 1].length;
          result.holes.push(holeIndex);
        }
      }
      return result;
    };
  }, {}] }, {}, [1]);
//# sourceMappingURL=scripts.js.map
