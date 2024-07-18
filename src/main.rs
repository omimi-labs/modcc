// mod ff;
// mod uni_poly;

use std::fmt::Debug;

use actix_web::{web, App, HttpResponse, HttpServer, Responder};
use serde::{Deserialize, Serialize};

struct LagrangeInput<T> {
    x_values: Vec<T>,
    y_values: Vec<T>,
    field: T,
}

#[actix_web::get("/lagrange_interpolation_over_ff/{user_id}/{user_name}/")]
async fn lagrange_interpolation_over_ff(
    path: web::Path<(String, String)>,
    // json: web::Json<LagrangeInput<T>>,
) -> impl Responder {
    let (user_id, user_name) = path.into_inner();
    println!("{}:{}", user_id, user_name);
    format!("{} {}", user_id, user_name)
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let port = 8080;
    println!("Server is running");

    HttpServer::new(move || App::new().service(lagrange_interpolation_over_ff))
        .bind(("127.0.0.1", port))?
        .workers(2)
        .run()
        .await
}
