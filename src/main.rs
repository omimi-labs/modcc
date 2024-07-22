mod ff;
mod interpolation;
mod uni_poly;

use actix_web::{get, post, web, App, HttpServer, Responder};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct LagrangeInput {
    x_values: Vec<usize>,
    y_values: Vec<usize>,
    field: usize,
}

#[get("/lf/")]
async fn lagrange_interpolation_over_ff(json: web::Json<LagrangeInput>) -> impl Responder {
    let x_values = &json.x_values;
    let y_values = &json.y_values;
    let field = json.field;
    // let poly = UniPoly::interpolate_xy(x_values, y_values);
    // println!("{:?}", poly.)
    "Lagrange View"
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
