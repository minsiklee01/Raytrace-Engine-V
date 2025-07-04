module main

import os
import math
import gfx

////////////////////////////////////////////////////////////////////////////////////////
// Comment out lines in array below to prevent re-rendering every scene.
// If you create a new scene file, add it to the list below.
//
// NOTE: **BEFORE** you submit your solution, uncomment all lines, so
//       your code will render all the scenes!

const scene_filenames = [
    'P02_00_sphere',
    'P02_01_sphere_ambient',
    'P02_02_sphere_room',
    'P02_03_quad',
    'P02_04_quad_room',
    'P02_05_ball_on_plane',
    'P02_06_balls_on_plane',
    'P02_07_reflections',
    'P02_08_antialiased',
    'P02_09_triangle',
    'P02_10_glass',
    'P02_11_ico_sphere'
]


////////////////////////////////////////////////////////////////////////////////////////
// module aliasing to make code a little easier to read
// ex: replacing `gfx.Scene` with just `Scene`

type Point     = gfx.Point
type Vector    = gfx.Vector
type Direction = gfx.Direction
type Normal    = gfx.Normal
type Ray       = gfx.Ray
type Color     = gfx.Color
type Image     = gfx.Image
type Shape     = gfx.Shape

type Intersection = gfx.Intersection
type Surface      = gfx.Surface
type Scene        = gfx.Scene
type Frame        = gfx.Frame


////////////////////////////////////////////////////////////////////////////////////////
// functions to implement


fn intersect_ray_surface(surface Surface, ray Ray) Intersection {

    if surface.shape == Shape.sphere {
        oc := ray.e - (surface.frame.o)
        a := ray.d.length_squared()
        b := 2.0 * (oc.dot(ray.d))
        c := oc.length_squared() - (surface.radius * surface.radius)
        d := b*b - 4.0*a*c

        if d < 0 { 
            return gfx.no_intersection
        }
    
        sqrt_d := math.sqrt(d)

        t1 := (-b - sqrt_d) / (2.0 * a)
        t2 := (-b + sqrt_d) / (2.0 * a)
        mut t := math.inf(1)

        if ray.valid_t(t1) && ray.valid_t(t2) {
            t = math.min(t1, t2)
        } else if ray.valid_t(t1) {
            t = t1
        } else if ray.valid_t(t2) {
            t = t2
        } else {
            return gfx.no_intersection
        }

        if !ray.valid_t(t) {
            return gfx.no_intersection
        }

        hit_point := ray.at(t)
        normal := (hit_point - surface.frame.o).normalize().as_direction()
        return Intersection{
            frame: gfx.Frame.new_oz(hit_point, normal)
            material: surface.material
            distance: t
        }
 
    } else if surface.shape == Shape.quad {
        n := surface.frame.z
        if ray.d.dot(n) == 0 {
            return gfx.no_intersection
        }

        t := (surface.frame.o - ray.e).dot(surface.frame.z) / ray.d.dot(n)

        if !ray.valid_t(t) {
            return gfx.no_intersection
        }

        hit_point := ray.at(t)
        local_point := surface.frame.w2l_point(hit_point)
        if math.abs(local_point.x) > surface.radius || math.abs(local_point.y) > surface.radius {
            return gfx.no_intersection
        }

        return Intersection{
            frame: gfx.Frame.new_oz(hit_point, surface.frame.z)
            material: surface.material
            distance: t
        }
    } else if surface.shape == Shape.triangle {
        ap := surface.frame.triangle_points.a - surface.frame.triangle_points.c
        bp := surface.frame.triangle_points.b - surface.frame.triangle_points.c
        ep := ray.e - surface.frame.triangle_points.c

        d_cross_bp := ray.d.cross(bp)
        denominator := d_cross_bp.dot(ap)

        t := (ep.cross(ap)).dot(bp) / denominator
        alpha := d_cross_bp.dot(ep) / denominator
        beta := ep.cross(ap).dot(ray.d) / denominator

        if ray.valid_t(t) && alpha >= 0 && beta >= 0 && (alpha + beta) <= 1 {
            hit_point := ray.at(t)
            normal := ap.cross(bp).normalize()
            return Intersection{
                frame: gfx.Frame.new_oz(hit_point, normal.as_direction())
                material: surface.material
                distance: t
            }
        } else {
            return gfx.no_intersection
        }
    }
    return gfx.no_intersection
}

// Determines if given ray intersects any surface in the scene.
// If ray does not intersect anything, null is returned.
// Otherwise, details of first intersection are returned as an `Intersection` struct.
fn intersect_ray_scene(scene Scene, ray Ray) Intersection {
    mut closest := gfx.no_intersection
        for surface in scene.surfaces {
            intersection := intersect_ray_surface(surface, ray)
            if intersection.hit() && intersection.is_closer(closest) {
                closest = intersection
            }
        }
        return closest
}

// Computes irradiance (as Color) from scene along ray
fn irradiance(scene Scene, ray Ray, depth int) Color {

    mut accum := gfx.black
    intersection := intersect_ray_scene(scene, ray)

    if intersection.miss() {
        return scene.background_color
    }

    // Ambient color
    accum += intersection.material.kd * scene.ambient_color

    for light in scene.lights {
        // light direction(l) and light response(L)
        l := (light.frame.o - intersection.o()).normalize()
        light_response := light.kl.lightness() / (light.frame.o - intersection.o()).length_squared()

        // Shadow ray to check visibility
        shadow_ray := Ray{
            e: intersection.o(),
            d: l.as_direction(),
            t_max: (light.frame.o - intersection.o()).magnitude()
        }
        // Visibility term (1 if not in shadow, 0 if in shadow)
        mut visibility := if intersect_ray_scene(scene, shadow_ray).miss() { 1.0 } else { 0.0 }

        // Blinn-Phong shading model
        kd := intersection.material.kd
        ks := intersection.material.ks

        view_dir := (ray.e - intersection.o()).normalize()
        h := (l + view_dir).normalize()
        brdf := ks.scale(math.pow(math.max(0.0, intersection.normal().dot(h)), intersection.material.n))
        c := (kd + brdf).scale(light_response * visibility * intersection.normal().dot(l))
        accum += c
    }

    // Reflection calculation
    if intersection.material.kr.lightness() > 0 {  
        reflection_ray := Ray{
            e: intersection.o(),
            d: ray.d.negate().reflect(intersection.normal())
            t_min: 10e-4
        }
        if depth - 1 >0 {
            reflection_color := irradiance(scene, reflection_ray, depth - 1)
            accum += reflection_color.mult(intersection.material.kr)
        }
    }

    // Refraction calculation
    if intersection.material.kt.is_black(){
        return accum
    }
    mut normal := intersection.normal()
    if normal.dot(ray.d) < 0.0 { 
        normal = normal.negate()
    }
    eta := 1.0 / intersection.material.kt.lightness() 
    c1 := normal.dot(ray.d.negate())
    c2_squared := eta * eta * (1 - c1 * c1)
    if c2_squared < 1.0 {
        c2 := math.sqrt(1 - c2_squared)
        // fresnel_result := fresnel(c1, c2 , intersection.material.kt.lightness())
        t := ray.d.negate().scale(-1.0 * eta) - normal.scale(eta*c1 + c2)
        refraction_ray := Ray{
            e: intersection.o(),
            d: t.as_direction()
        }
        if depth - 1 >0 {
            refraction_color := irradiance(scene, refraction_ray, depth - 1)
            accum += refraction_color.scale(1.5)
            // .scale(fresnel_result.transmittance)
        }
    }
    
    return accum
}
pub struct FresnelResult {
    reflectance f64
    transmittance f64
}
/*
fn fresnel(cos_theta1 f64, cos_theta2, n_to f64) FresnelResult {
    fr1 := math.pow(( n_to * cos_theta1 - 1.0 * cos_theta2 ) / ( n_to * cos_theta1 + 1.0 * cos_theta2), 2)
    // fr1 := math.pow(( 1.0 * cos_theta1 - n_to * cos_theta2 ) / ( 1.0 * cos_theta1 + n_to * cos_theta2), 2)
    fr2 := math.pow(( 1.0 * cos_theta2 - n_to * cos_theta1 ) / ( 1.0 * cos_theta2 + n_to * cos_theta1), 2)
    // fr2 := math.pow(( n_to * cos_theta2 - 1.0 * cos_theta1 ) / ( n_to * cos_theta2 + 1.0 * cos_theta1), 2)
    fr := (fr1 + fr2) * 0.5
    ft := 1 - fr
    return FresnelResult{reflectance: fr, transmittance: ft}
}
*/
fn raytrace(scene Scene) Image {
    // convenience vars
    camframe   := scene.camera.frame
    sensor     := scene.camera.sensor
    resolution := sensor.resolution
    size       := sensor.size

    mut image := gfx.Image.new(resolution)
    /*
    for y in 0 .. resolution.height {
        for x in 0 .. resolution.width {
            // convert (x,y) to normalized coords (u,v)
            u := (f64(x) + 0.5) / f64(resolution.width)
            v := (f64(resolution.height) - (f64(y) + 0.5)) / f64(resolution.height)

            // compute image plane offsets
            fx := camframe.x.scale((u - 0.5) * size.width)
            fy := camframe.y.scale((v - 0.5) * size.height)
            fz := camframe.z.scale(-sensor.distance)

            // compute point at center of pixel in image plane
            q := camframe.o.add(fx).add(fy).add(fz)
            ray := camframe.o.ray_through(q)
            color := irradiance(scene, ray, 30)
            image.set_xy(x, y, color) 
        }
    }
    */
    mut color := gfx.black
    s := scene.camera.sensor.samples
    for y in 0 .. resolution.height {
        for x in 0 .. resolution.width {
            color = gfx.black
            for i in 0 .. s {
                for j in 0 .. s{
                    u := (f64(x) + (f64(i) + 0.5) / f64(s)) / f64(resolution.width)
                    v := (f64(resolution.height) - f64(y) - (f64(j) + 0.5) / f64(s)) / f64(resolution.height)

                    // compute image plane offsets
                    fx := camframe.x.scale((u - 0.5) * size.width)
                    fy := camframe.y.scale((v - 0.5) * size.height)
                    fz := camframe.z.scale(-sensor.distance)

                    // compute point at center of pixel in image plane
                    q := camframe.o.add(fx).add(fy).add(fz)
                    ray := camframe.o.ray_through(q)
                    color.add_in(irradiance(scene, ray, 5)) 
                }
            }
            image.set_xy(x, y, color.scale(1.0/f64(s*s))) 
        }
    }
    
    return image
    
}

fn main() {
    // Make sure images folder exists, because this is where all generated images will be saved
    if !os.exists('output') {
        os.mkdir('output') or { panic(err) }
    }

    for filename in scene_filenames {
        println('Rendering ${filename}...')
        scene := gfx.scene_from_file('scenes/${filename}.json')!
        image := raytrace(scene)
        image.save_png('output/${filename}.png')
    }

    println('Done!')
}
