mod ff;
mod lagrange_interpolation;
mod multilinear_interpolation;
mod multilinear_poly;
mod univariate_poly;

use univariate_poly::LagrangeInterpolationSteps;

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

#[derive(Debug, Deserialize)]
struct MultilinearOverBooleanHypercubeRequest {
    y_values: Vec<i128>,
    field: usize,
}

#[derive(Debug, Serialize)]
struct MultilinearOverBooleanHypercubeResponse {
    coefficients: Vec<(usize, usize)>,
}

#[post("/lagrange_interpolation/")]
async fn lagrange_interpolation_over_ff(
    json: web::Json<LagrangeInterpolationRequest>,
) -> impl Responder {
    let x_values = &json.x_values;
    let y_values = &json.y_values;
    let field = json.field;
    let (coefficients, steps) =
        lagrange_interpolation::lagrange_interpolate(x_values, y_values, field as u128);
    let response = LagrangeInterpolationResponse {
        coefficients,
        steps,
    };
    HttpResponse::Ok().json(response)
}

#[post("/mf/")]
async fn multilinear_interpolation_over_boolean_hypercube(
    json: web::Json<MultilinearOverBooleanHypercubeRequest>,
) -> impl Responder {
    let y_values = &json.y_values;
    let field = json.field;
    let coefficients = multilinear_interpolation::multilinear_interpolate_over_boolean_hypercube(
        y_values,
        field as u128,
    );
    let response = MultilinearOverBooleanHypercubeResponse { coefficients };
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
