mod ff;
mod interpolation;
mod uni_poly;

use uni_poly::LagrangeInterpolationSteps;

use actix_cors::Cors;
use actix_web::{http, post, web, App, HttpResponse, HttpServer, Responder};
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize)]
struct LagrangeInterpolationRequest {
    x_values: Vec<i128>,
    y_values: Vec<i128>,
    field: usize,
}

#[derive(Debug, Serialize)]
struct LagrangeInterpolationResponse {
    coefficients: Vec<u128>,
    steps: LagrangeInterpolationSteps,
}

#[post("/lf/")]
async fn lagrange_interpolation_over_ff(
    json: web::Json<LagrangeInterpolationRequest>,
) -> impl Responder {
    let x_values = &json.x_values;
    let y_values = &json.y_values;
    let field = json.field;
    let (coefficients, steps) =
        interpolation::lagrange_interpolate(x_values, y_values, field as u128);
    let response = LagrangeInterpolationResponse {
        coefficients,
        steps,
    };
    println!("{:?}", response);
    HttpResponse::Ok().json(response)
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let port = 8080;
    println!("Server is running");

    HttpServer::new(move || {
        App::new()
            .wrap(
                Cors::default()
                    .allowed_origin("http://localhost:3000") // Specify your Next.js app's origin
                    .allowed_methods(vec!["GET", "POST"])
                    .allowed_headers(vec![http::header::AUTHORIZATION, http::header::ACCEPT])
                    .allowed_header(http::header::CONTENT_TYPE)
                    .max_age(3600),
            )
            .service(lagrange_interpolation_over_ff)
    })
    .bind(("127.0.0.1", port))?
    .workers(2)
    .run()
    .await
}
