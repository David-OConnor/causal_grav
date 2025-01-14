//! https://github.com/BaGreal2/gravitation-particles/tree/main/src

use std::fs;
use lin_alg::f64::Vec3;
use crate::Body;

#[derive(Clone)]
pub struct OctTree {
    bounds: Cube,
    children: Option<[Box<OctTree>; 8]>,
    body: Option<Body>,
    mass: f64,
    m_center_pos: Vec3,
}

impl OctTree {
    pub fn new(bounds: Cube) -> Self {
        let copy_bounds = bounds.clone();
        Self {
            bounds,
            children: None,
            body: None,
            mass: 0.0,
            m_center_pos: Vec3::new(
                                        copy_bounds.top_left_fwd_pos.x + copy_bounds.w / 2.0,
                                        copy_bounds.top_left_fwd_pos.y + copy_bounds.h / 2.0,
                                        copy_bounds.top_left_fwd_pos.z + copy_bounds.d / 2.0,
            )
        }
    }

    fn is_divided(&self) -> bool {
        !self.children.is_none()
    }

    fn subdivide(&mut self) {
        let (x, y, z) = (self.bounds.top_left_fwd_pos.x, self.bounds.top_left_fwd_pos.y, self.bounds.top_left_fwd_pos.z);
        let (w, h, d) = (self.bounds.w, self.bounds.h, self.bounds.d);

        let top_left_f = Cube::new(Vec3::new(x, y, z), w / 2.0, h / 2.0, d/2.);
        let top_right_f = Cube::new(Vec3::new(x + w / 2.0, y, z), w / 2.0, h / 2.0, d/2.);
        let bottom_left_f = Cube::new(Vec3::new(x, y + h / 2.0, z), w / 2.0, h / 2.0, d/2.);
        let bottom_right_f = Cube::new(Vec3::new(x + w / 2.0, y + h / 2.0, z), w / 2.0, h / 2.0, d/2.);

        let top_left_a = Cube::new(Vec3::new(x, y, z), w / 2.0, h / 2.0, d/2.);
        let top_right_a = Cube::new(Vec3::new(x + w / 2.0, y, z + d/2.), w / 2.0, h / 2.0, d/2.);
        let bottom_left_a = Cube::new(Vec3::new(x, y + h / 2.0, z + d/2.), w / 2.0, h / 2.0, d/2.);
        let bottom_right_a = Cube::new(Vec3::new(x + w / 2.0, y + h / 2.0, z + d/2.), w / 2.0, h / 2.0, d/2.);
        
        self.children = Some([
            Box::new(OctTree::new(top_left_f)),
            Box::new(OctTree::new(top_right_f)),
            Box::new(OctTree::new(bottom_left_f)),
            Box::new(OctTree::new(bottom_right_f)),
            Box::new(OctTree::new(top_left_a)),
            Box::new(OctTree::new(top_right_a)),
            Box::new(OctTree::new(bottom_left_a)),
            Box::new(OctTree::new(bottom_right_a)),
        ]);
    }

    pub fn insert(&mut self, particle: &Body) {
        if !self.bounds.contains(particle) {
            return;
        }

        if self.body.is_none() {
            self.body = Some(*particle);
        } else {
            if !self.is_divided() {
                self.subdivide();
            }
            // self.children
            //     .as_mut()
            //     .unwrap()
            //     .par_iter_mut()
            //     .for_each(|leaf| {
            //         leaf.as_mut().insert(particle);
            //     });
            for leaf in self.children.as_mut().unwrap().as_mut() {
                leaf.as_mut().insert(particle);
            }
            self.update_mass();
        }
    }

    pub fn calculate_force(&mut self, particle: &mut Particle) {
        if !self.is_divided() {
            if let Some(existent_particle) = &self.body {
                if existent_particle.index != particle.index {
                    let attraction_force =
                        particle.get_attraction_force(&self.body.as_ref().unwrap());
                    // let attraction_force = self.particle.unwrap().get_attraction_force(particle);
                    particle.net_force += attraction_force;
                }
            }
            return;
        }

        let ratio = self.bounds.w / particle.get_distance_to(&self.m_center_pos);
        if ratio < 0.5 {
            let attraction_force = particle.get_attraction_force(&Body::new(
                self.m_center_pos,
                Vec3::new(0.0, 0.0),
                self.mass,
                1.0,
                1000000,
            ));
            particle.net_force += attraction_force;
            return;
        }

        for leaf in self.children.as_mut().unwrap() {
            leaf.calculate_force(particle);
        }
    }

    fn update_mass(&mut self) {
        if !self.is_divided() {
            if self.body.is_none() {
                return;
            }
            self.mass = self.body.as_ref().unwrap().mass;
            self.m_center_pos = self.body.as_ref().unwrap().pos;
            return;
        }
        let mut mass_sum = 0.0;
        let mut center_x= 0.0;
        let mut center_y = 0.0;
        let mut center_z = 0.0;

        for leaf in self.children.as_mut().unwrap() {
            leaf.as_mut().update_mass();
            mass_sum += leaf.as_ref().mass;

            center_x += leaf.as_ref().m_center_pos.x * leaf.as_ref().mass;
            center_y += leaf.as_ref().m_center_pos.y * leaf.as_ref().mass;
            center_z += leaf.as_ref().m_center_pos.z * leaf.as_ref().mass;
        }
        self.mass = mass_sum;

        center_x /= mass_sum;
        center_y /= mass_sum;
        center_z /= mass_sum;

        self.m_center_pos = Vec3::new(center_x, center_y, center_z);
    }

    pub fn show(
        &self,
        canvas: &mut Canvas,
        ctx: &mut Context,
        offset: Vector2<f64>,
        zoom: f64,
        particles_to_draw: &Vec<Body>,
        max_vel: f64,
        min_vel: f64,
        show_bounds: bool,
    ) {
        if show_bounds {
            self.bounds.show(
                canvas,
                ctx,
                offset,
                zoom,
                &mut Color::from_rgb(255, 255, 255),
            );
        }
        match self.children.as_ref() {
            Some(children) => {
                for leaf in children.as_ref() {
                    leaf.show(
                        canvas,
                        ctx,
                        offset,
                        zoom,
                        particles_to_draw,
                        max_vel,
                        min_vel,
                        show_bounds,
                    );
                }
            }
            None => {}
        }

        match &self.body {
            Some(existent_particle) => {
                if particles_to_draw.contains(&existent_particle) {
                    existent_particle.show(canvas, ctx, offset, zoom, max_vel, min_vel);
                }
            }
            None => {}
        }
    }

    pub fn query(&self, rect: &Cube) -> Vec<Body> {
        let mut results = Vec::new();
        if !self.bounds.intersects(rect) {
            return results;
        }

        if let Some(particle) = &self.body {
            if rect.contains(particle) {
                results.push(*particle);
            }
        }
        if self.is_divided() {
            for leaf in self.children.as_ref().unwrap().iter() {
                results.extend(leaf.query(rect));
            }
        }

        results
    }
}

fn random_in_circle(radius: f64, padding: f64, center: Vector2<f64>) -> Vector2<f64> {
    let mut rng = rand::thread_rng();
    let angle = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
    let distance = rng.gen_range(padding..radius);

    Vec3::new(distance * angle.cos(), distance * angle.sin()) + center
}

pub fn spawn_circle(
    particles: &mut Vec<Body>,
    center: Vector2<f64>,
    radius: f64,
    particle_mass: f64,
    particles_amount: i32,
) {
    for i in 0..particles_amount {
        let pos = random_in_circle(radius, 0.0, center);
        let new_particle =
            Particle::new(pos, Vector2::default(), particle_mass, 0.00001, i as usize);
        particles.push(new_particle);
    }
}

pub fn create_galaxy(
    particles: &mut Vec<Body>,
    center: Vec3,
    initial_vel: Vector2<f64>,
    radius: f64,
    sun_mass: f64,
    particle_mass: f64,
    particles_amount: i32,
) {
    for i in 0..particles_amount {
        let pos = random_in_circle(radius, 2.0, center);
        let distance_to_center = pos.metric_distance(&center);
        let orbital_vel = ((G * sun_mass) / distance_to_center).sqrt();
        let dir = Vec3::new(pos.y - center.y, center.x - pos.x).to_normalized();
        let new_particle =
            Particle::new(pos, dir * orbital_vel, particle_mass, 0.00001, i as usize);
        particles.push(new_particle);
    }

    let sun = Particle::new(
        center,
        initial_vel,
        sun_mass,
        1.5,
        particles_amount as usize,
    );
    particles.push(sun);
}

pub fn create_quadtree(particles: &Vec<Body>) -> OctTree {
    let mut qt = OctTree::new(Cube::new(
        Vec3::new(0., 0., 0.),
        WORLD_WIDTH,
        WORLD_HEIGHT,
    ));
    for i in 0..particles.len() {
        qt.insert(&Bodys[i]);
    }
    qt
}

pub fn calculate_new_position(particle: &mut Particle, qt: &mut OctTree) {
    particle.net_force = Vec3::new(0.0, 0.0);
    qt.calculate_force(particle);
    // println!("{:?}", borrowed.net_force);

    let acceleration = particle.net_force / particle.mass;
    particle.vel += acceleration;
    let velocity = particle.vel;
    particle.pos += velocity;
}

pub fn world_to_screen_coords(
    world_coords: Vector2<f64>,
    origin: &Vector2<f64>,
    zoom: f64,
) -> Vector2<f64> {
    (origin + world_coords) * zoom
}
pub fn screen_to_world_coords(
    screen_coords: Vector2<f64>,
    origin: &Vector2<f64>,
    zoom: f64,
) -> Vector2<f64> {
    screen_coords / zoom - origin
}

pub fn rename_images(ctx: &Context) {
    let data_dir = ctx.fs.user_data_dir();
    for file in fs::read_dir(data_dir.join("image-cache/")).unwrap() {
        let full_path: PathBuf = file.as_ref().unwrap().path();
        let full_path_string: String = full_path.to_string_lossy().to_string();
        let full_name = String::from(
            file.as_ref()
                .unwrap()
                .path()
                .file_name()
                .unwrap()
                .to_string_lossy(),
        );
        let name = full_name[0..full_name.len() - 4].to_owned();
        if name.chars().nth(0).unwrap() == '.' {
            continue;
        }
        let prefix_amount = 6 - name.len();
        let repeated_string: String = std::iter::repeat("0").take(prefix_amount).collect();
        let path = &full_path_string[0..full_path_string.len() - full_name.len()];
        let old_path = String::from(path) + &full_name;
        let new_path = String::from(path) + &repeated_string + &full_name;
        let _ = fs::rename(old_path, new_path);
    }
}

pub fn convert_to_video(ctx: &Context) {
    let data_dir = ctx.fs.user_data_dir().to_string_lossy().to_string();
    let local: DateTime<Local> = Local::now();
    let formatted_date_time = local.format("%Y-%m-%d_%H.%M.%S").to_string();
    let mut cmd = Command::new("ffmpeg")
        .args([
            "-framerate",
            "60",
            "-pattern_type",
            "glob",
            "-i",
            format!("{data_dir}/image-cache/*.png").as_str(),
            "-vf",
            "eq=saturation=2",
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            format!("results/{formatted_date_time}.mp4").as_str(),
        ])
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    {
        let stdout = cmd.stdout.as_mut().unwrap();
        let stdout_reader = BufReader::new(stdout);
        let stdout_lines = stdout_reader.lines();

        for line in stdout_lines {
            println!("Read: {:?}", line);
        }
    }

    cmd.wait().unwrap();
}

pub fn clean_cache_images(ctx: &Context) {
    let data_dir = ctx.fs.user_data_dir();
    for file in fs::read_dir(data_dir.join("image-cache/")).unwrap() {
        let full_path: PathBuf = file.as_ref().unwrap().path();
        let _ = fs::remove_file(full_path);
    }
}

pub fn move_on_mouse(ctx: &mut Context, origin: &mut Vector2<f64>, zoom: f64) {
    const DESIRED_FPS: u32 = 60;

    while ctx.time.check_update_time(DESIRED_FPS) {
        let mouse_position = ctx.mouse.position();

        if mouse_position.x < LOWER_BOUND.x {
            origin.x += 5.0;
        } else if mouse_position.x > UPPER_BOUND.x {
            origin.x -= 5.0;
        }
        if mouse_position.y < LOWER_BOUND.y {
            origin.y += 5.0;
        } else if mouse_position.y > UPPER_BOUND.y {
            origin.y -= 5.0;
        }

        if -origin.x < 0.0 {
            origin.x = 0.0;
        } else if -origin.x + WIDTH / zoom > WORLD_WIDTH {
            origin.x = WIDTH / zoom - WORLD_WIDTH;
        }
        if -origin.y < 0.0 {
            origin.y = 0.0;
        } else if -origin.y + HEIGHT / zoom > WORLD_HEIGHT {
            origin.y = HEIGHT / zoom - WORLD_HEIGHT;
        }
    }
}

pub fn zoom_world(ctx: &Context, origin: &mut Vector2<f64>, zoom: &mut f64, y_diff: f64) {
    let mouse_x = ctx.mouse.position().x;
    let mouse_y = ctx.mouse.position().y;

    let mut mouse_world = screen_to_world_coords(Vec3::new(mouse_x, mouse_y), origin, *zoom);
    if y_diff > 0.0 {
        *zoom *= 1.1;
    } else if y_diff < 0.0 {
        *zoom /= 1.1;
    }
    if *zoom < MAX_ZOOM {
        *zoom = MAX_ZOOM;
    }
    mouse_world = world_to_screen_coords(mouse_world, origin, *zoom);

    origin.x += (mouse_x - mouse_world.x) / *zoom;
    origin.y += (mouse_y - mouse_world.y) / *zoom;
}

pub fn save_screen(ctx: &Context, screen: &mut ScreenImage, frame_count: u32) {
    let output_name = String::from("/image-cache/") + frame_count.to_string().as_str() + ".png";
    match screen
        .image(ctx)
        .encode(ctx, ImageEncodingFormat::Png, output_name)
    {
        Err(saving_err) => eprintln!("{}", saving_err),
        _ => {}
    }
}


#[derive(Clone)]
pub struct Cube {
    pub top_left_fwd_pos: Vec3,
    pub w: f64,
    pub h: f64,
    pub d: f64,
}

impl Cube {
    pub fn new(top_left_pos: Vector2<f64>, w: f64, h: f64, d: f64) -> Self {
        Self { top_left_fwd_pos: top_left_pos, w, h, d }
    }

    pub fn contains(&self, particle: &Body) -> bool {
        self.top_left_fwd_pos.x <= particle.pos.x
            && self.top_left_fwd_pos.y <= particle.pos.y
            && self.top_left_fwd_pos.x + self.w > particle.pos.x
            && self.top_left_fwd_pos.y + self.h > particle.pos.y
    }

    pub fn intersects(&self, rect: &Self) -> bool {
        let up = rect.top_left_fwd_pos.y + rect.h < self.top_left_fwd_pos.y;
        let down = rect.top_left_fwd_pos.y > self.top_left_fwd_pos.y + self.h;
        let left = rect.top_left_fwd_pos.x + rect.w < self.top_left_fwd_pos.x;
        let right = rect.top_left_fwd_pos.x > self.top_left_fwd_pos.x + self.w;
        !(up || down || left || right)
    }

    pub fn show(
        &self,
        canvas: &mut Canvas,
        ctx: &mut Context,
        offset: Vector2<f64>,
        zoom: f64,
        color: &mut Color,
    ) {
        color.a = 0.3;
        let rect = graphics::Rect {
            x: world_to_screen_coords(self.top_left_fwd_pos, &offset, zoom).x,
            y: world_to_screen_coords(self.top_left_fwd_pos, &offset, zoom).y,
            w: self.w * zoom,
            h: self.h * zoom,
        };
        let rect_mesh = graphics::Mesh::new_rectangle(
            ctx,
            graphics::DrawMode::Stroke(graphics::StrokeOptions::DEFAULT),
            rect,
            *color,
        )
            .unwrap();

        canvas.draw(&rect_mesh, graphics::DrawParam::default());
    }
}