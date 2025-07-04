module gfx

import os
import json


///////////////////////////////////////////////////////////
// structs used to define a scene
//
// NOTE: if any post-load processing needs to happen on
//       structs after loading from JSON file, add
//       processing code to `update_after_load` below!
//
//       ex: Frames can have Look At parameters (eye, target, up)
//           specified instead of usual Frame params (o, x, y, z).

pub struct Scene {
pub:
    camera           Camera
    background_color Color       = Color{ 0.2, 0.2, 0.2 }
    ambient_color    Color       = Color{ 0.2, 0.2, 0.2 }
    lights           [] Light    = [ Light{} ]
    surfaces         [] Surface  = [ Surface{} ]
}

pub struct Camera {
pub:
    sensor Sensor = Sensor{}
    frame  Frame  = Frame{
        o: Point{ 0.0, 0.0, 1.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
}

// Sensor is also called Image Plane
pub struct Sensor {
pub:
    size       Size2  = Size2{  1.0, 1.0 }  // "physical" size of sensor
    resolution Size2i = Size2i{ 512, 512 }  // resolution of sensor (pixels count)
    distance   f64    = 1.0                 // distance between camera focal point and sensor
    samples    int    = 1                   // samples per pixel
}

pub struct Light {
pub:
    kl     Color = white                    // brightness and color of light
    frame  Frame = Frame{
        o: Point{ 0.0, 0.0, 5.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
}

pub enum Shape {
    sphere
    quad
    triangle
}


pub struct Surface {
pub:
    shape    Shape = Shape.sphere           // shape of surface: sphere or quad
    radius   f64   = 1.0                    // size of surface
    frame    Frame
    material Material
}

pub struct Material {
pub:
    kd Color = white                         // diffuse reflection
    ks Color = black                         // blinn-phong specular reflection 
    n  f64   = 10.0                          // blinn-phong specular exponent
    kr Color = black                         // mirror reflection
    kt Color = black // transmission coefficient
}

pub struct TrianglePoints {
pub:
    a Point =  Point{ -2.5, -1.25, 0.0}
    b Point = Point{ 2.5, -1.25, 0.0}
    c Point = Point{ 0, 1.25, 0.0}
}
// function to update scene after loading from JSON
fn (scene Scene) update_after_load() Scene {
    // NOTE: attempted a more generic version using reflection (see bottom)
    return Scene{
        ...scene                                            // copy everything from scene with some overrides (below)
        camera: Camera{
            ...scene.camera                                 // copy everything from camera with override
            frame: scene.camera.frame.as_frame()            // update frame for camera if look at is specified
        }
        lights: scene.lights.map(fn (light Light) Light {   // for each light...
            return Light{
                ...light                                    // copy everything from light with override
                frame: light.frame.as_frame()               // update frame for light if look at is specified
            }
        })
        surfaces: scene.surfaces.map(fn (surface Surface) Surface {
            return Surface{
                ...surface                                  // copy everything from surface with override
                frame: surface.frame.as_frame()             // update frame for surface if look at is specified
            }
        })
    }
}



///////////////////////////////////////////////////////////
// scene importing and exporting functions

pub fn scene_from_file(path string) !Scene {
    // load and decode scene from JSON file
    data := os.read_file(path)!
    scene := json.decode(Scene, data)!
    return scene.update_after_load()
}

pub fn (scene Scene) to_json() string {
    return json.encode(scene)
}



///////////////////////////////////////////////////////////
// convenience getters

pub fn (light Light) o() Point {
    return light.frame.o
}

pub fn (surface Surface) o() Point {
    return surface.frame.o
}



///////////////////////////////////////////////////////////
// attempting a more generic updated function that uses reflection

// // commented out following code due to the wrong struct data is seen on the lines marked with <---
// type SceneType = Scene | Camera | Sensor | Light | Surface | Material
// fn update_frame_from_load(f Frame) Frame {
//     println('UPDATING FRAME FROM LOAD')
//     println(f)                                          // <---  this has INCORRECT data!!
//     if f.eye != f.target && f.up != direction_no {
//         // look at was specified!
//         ret := frame_lookat(f.eye, f.target, f.up)
//         println(ret)
//         return ret
//     }
//     return f
// }
// fn update_array_from_load[T](a []T) []T {
//     mut result := []T{}
//     for i in a {
//         result << update_after_load(i)
//     }
//     return result
// }
// fn update_after_load[T](o T) T {
//     mut result := T{}
//     $for field in T.fields {
//         $if field.typ is $array {
//             result.$(field.name) = update_array_from_load(o.$(field.name))
//         } $else $if field.typ is $enum {
//             result.$(field.name) = o.$(field.name)
//         } $else $if field.typ is Frame {
//             println(o.$(field.name))                                        // <---  this has INCORRECT data!!
//             result.$(field.name) = update_frame_from_load(o.$(field.name))
//         } $else $if field.typ is SceneType {
//             result.$(field.name) = update_after_load(o.$(field.name))
//         } $else {
//             result.$(field.name) = o.$(field.name)
//         }
//     }
//     return result
// }





