// Provides a convenient raytracing interface.

// ### new GL.HitTest([t, hit, normal])
//
// This is the object used to return hit test results. If there are no
// arguments, the constructed argument represents a hit infinitely far
// away.
function HitTest(t, hit, normal) {
  this.t = arguments.length ? t : Number.MAX_VALUE;
  this.hit = hit;
  this.normal = normal;
}

// ### .mergeWith(other)
//
// Changes this object to be the closer of the two hit test results.
HitTest.prototype = {
  mergeWith: function(other) {
    if (other.t > 0 && other.t < this.t) {
      this.t = other.t;
      this.hit = other.hit;
      this.normal = other.normal;
    }
  }
};

// ### new GL.Raytracer()
//
// This will read the current modelview matrix, projection matrix, and viewport,
// reconstruct the eye position, and store enough information to later generate
// per-pixel rays using `getRayForPixel()`.
//
// Example usage:
//
//     var tracer = new GL.Raytracer();
//     var ray = tracer.getRayForPixel(
//       gl.canvas.width / 2,
//       gl.canvas.height / 2);
//     var result = GL.Raytracer.hitTestSphere(
//       tracer.eye, ray, new GL.Vector(0, 0, 0), 1);
function Raytracer() {
  var v = gl.getParameter(gl.VIEWPORT);
  var m = gl.modelviewMatrix.m;

  var axisX = new Vector(m[0], m[4], m[8]);
  var axisY = new Vector(m[1], m[5], m[9]);
  var axisZ = new Vector(m[2], m[6], m[10]);
  var offset = new Vector(m[3], m[7], m[11]);
  this.eye = new Vector(-offset.dot(axisX), -offset.dot(axisY), -offset.dot(axisZ));

  var minX = v[0], maxX = minX + v[2];
  var minY = v[1], maxY = minY + v[3];
  this.ray00 = gl.unProject(minX, minY, 1).subtract(this.eye);
  this.ray10 = gl.unProject(maxX, minY, 1).subtract(this.eye);
  this.ray01 = gl.unProject(minX, maxY, 1).subtract(this.eye);
  this.ray11 = gl.unProject(maxX, maxY, 1).subtract(this.eye);
  this.viewport = v;
}

Raytracer.prototype = {
  // ### .getRayForPixel(x, y)
  //
  // Returns the ray originating from the camera and traveling through the pixel `x, y`.
  getRayForPixel: function(x, y) {
    x = (x - this.viewport[0]) / this.viewport[2];
    y = 1 - (y - this.viewport[1]) / this.viewport[3];
    var ray0 = Vector.lerp(this.ray00, this.ray10, x);
    var ray1 = Vector.lerp(this.ray01, this.ray11, x);
    return Vector.lerp(ray0, ray1, y).unit();
  }
};

// ### GL.Raytracer.hitTestBox(origin, ray, pos, orientationMatrix, size)
//
// Traces the ray starting from `origin` along `ray` against the object bouding box
// whose size is `size`, orientation is `orientationMatrix` and position is `pos`. Returns a `HitTest` with the
// information or `null` for no intersection.
Raytracer.hitTestOBB = function(origin, ray, pos, orientationMatrix, size) {
    const p = pos.subtract(origin);
    const m = orientationMatrix.m;
    const x = new Vector(m[0], m[1], m[2]);
    const y = new Vector(m[4], m[5], m[6]);
    const z = new Vector(m[8], m[9], m[10]);

    const s = [size.x, size.y, size.z];
    const f = [x.dot(ray), y.dot(ray), z.dot(ray)];
    const e = [x.dot(p), y.dot(p), z.dot(p)];
    const t = [0, 0, 0, 0, 0, 0];
    for(let i = 0; i < 3; i ++){
      if(Math.abs(f[i]) < 1.0e-6){
        return null;
      }
      t[i * 2 + 0] = (e[i] + s[i]) / f[i];
      t[i * 2 + 1] = (e[i] - s[i]) / f[i];
    }

    let tmin = Math.max(Math.max(Math.min(t[0], t[1]), Math.min(t[2], t[3])), Math.min(t[4], t[5]));
    let tmax = Math.min(Math.min(Math.max(t[0], t[1]), Math.max(t[2], t[3])), Math.max(t[4], t[5]));

    if(tmax < 0 || tmin > tmax)
      return null;

    let time = tmin;
    if(tmin < 0)
      time = tmax;

    return new HitTest(time, origin.add(ray.multiply(time)), new GL.Vector(1, 0, 0));
};

// ### GL.Raytracer.hitTestBox(origin, ray, min, max)
//
// Traces the ray starting from `origin` along `ray` against the axis-aligned box
// whose coordinates extend from `min` to `max`. Returns a `HitTest` with the
// information or `null` for no intersection.
//
// This implementation uses the [slab intersection method](http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm).
Raytracer.hitTestBox = function(origin, ray, min, max) {
  var tMin = min.subtract(origin).divide(ray);
  var tMax = max.subtract(origin).divide(ray);
  var t1 = Vector.min(tMin, tMax);
  var t2 = Vector.max(tMin, tMax);
  var tNear = t1.max();
  var tFar = t2.min();

  if (tNear > 0 && tNear < tFar) {
    var epsilon = 1.0e-6, hit = origin.add(ray.multiply(tNear));
    min = min.add(epsilon);
    max = max.subtract(epsilon);
    return new HitTest(tNear, hit, new Vector(
      (hit.x > max.x) - (hit.x < min.x),
      (hit.y > max.y) - (hit.y < min.y),
      (hit.z > max.z) - (hit.z < min.z)
    ));
  }

  return null;
};

// ### GL.Raytracer.hitTestSphere(origin, ray, center, radius)
//
// Traces the ray starting from `origin` along `ray` against the sphere defined
// by `center` and `radius`. Returns a `HitTest` with the information or `null`
// for no intersection.
Raytracer.hitTestSphere = function(origin, ray, center, radius) {
  var offset = origin.subtract(center);
  var a = ray.dot(ray);
  var b = 2 * ray.dot(offset);
  var c = offset.dot(offset) - radius * radius;
  var discriminant = b * b - 4 * a * c;

  if (discriminant > 0) {
    var t = (-b - Math.sqrt(discriminant)) / (2 * a), hit = origin.add(ray.multiply(t));
    return new HitTest(t, hit, hit.subtract(center).divide(radius));
  }

  return null;
};

// ### GL.Raytracer.hitTestTriangle(origin, ray, a, b, c)
//
// Traces the ray starting from `origin` along `ray` against the triangle defined
// by the points `a`, `b`, and `c`. Returns a `HitTest` with the information or
// `null` for no intersection.
Raytracer.hitTestTriangle = function(origin, ray, a, b, c) {
  var ab = b.subtract(a);
  var ac = c.subtract(a);
  var normal = ab.cross(ac).unit();
  var t = normal.dot(a.subtract(origin)) / normal.dot(ray);

  if (t > 0) {
    var hit = origin.add(ray.multiply(t));
    var toHit = hit.subtract(a);
    var dot00 = ac.dot(ac);
    var dot01 = ac.dot(ab);
    var dot02 = ac.dot(toHit);
    var dot11 = ab.dot(ab);
    var dot12 = ab.dot(toHit);
    var divide = dot00 * dot11 - dot01 * dot01;
    var u = (dot11 * dot02 - dot01 * dot12) / divide;
    var v = (dot00 * dot12 - dot01 * dot02) / divide;
    if (u >= 0 && v >= 0 && u + v <= 1) return new HitTest(t, hit, normal);
  }

  return null;
};
