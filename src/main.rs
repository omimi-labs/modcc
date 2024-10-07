mod ff;
mod multilinear_poly;
mod multilinear_view_helpers;
mod multivariate_poly;
mod multivariate_view_helpers;
mod univariate_poly;
mod univariate_view_helpers;

use multilinear_poly::MultilinearLagrangeInterpolationSteps;
use multivariate_view_helpers::multivariate_interpolate_over_finite_field;
use univariate_poly::LagrangeInterpolationSteps;

use actix_cors::Cors;
use actix_web::{http, post, web, App, HttpResponse, HttpServer, Responder};
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize)]
struct EvaluateUnivariatePolyRequest {
    evaluation_point: isize,
    poly_string: String,
    field: usize,
}

#[derive(Debug, Serialize)]
struct EvaluateUnivariatePolyResponse {
    evaluation: usize,
}

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
    steps: MultilinearLagrangeInterpolationSteps,
}

#[derive(Debug, Deserialize)]
struct MultivariateInterpolationRequest {
    num_of_vars: usize,
    evaluation_points: Vec<Vec<i128>>,
    y_values: Vec<i128>,
    field: usize,
}

#[derive(Debug, Serialize)]
struct MultivariateInterpolationResponse {
    terms: Vec<(usize, Vec<(usize, usize)>)>,
}

#[derive(Serialize)]
struct BadRequestResponse {
    detail: String,
}

#[post("/evaluate_univariate_poly/")]
async fn evaluate_univariate_poly(
    json: web::Json<EvaluateUnivariatePolyRequest>,
) -> impl Responder {
    let poly_string = &json.poly_string;
    let field = json.field;
    let evaluation_point = json.evaluation_point;
    let evaluation = univariate_view_helpers::evaluate(evaluation_point, &poly_string, field);
    let response = EvaluateUnivariatePolyResponse { evaluation };
    HttpResponse::Ok().json(response)
}

#[post("/lagrange_interpolation/")]
async fn lagrange_interpolation_over_ff(
    json: web::Json<LagrangeInterpolationRequest>,
) -> impl Responder {
    let x_values = &json.x_values;
    let y_values = &json.y_values;
    let field = json.field;
    let (coefficients, steps) =
        univariate_view_helpers::lagrange_interpolate(x_values, y_values, field as u128);
    let response = LagrangeInterpolationResponse {
        coefficients,
        steps,
    };
    HttpResponse::Ok().json(response)
}

#[post("/multilinear_interpolation_over_boolean_hypercube/")]
async fn multilinear_interpolation_over_boolean_hypercube(
    json: web::Json<MultilinearOverBooleanHypercubeRequest>,
) -> impl Responder {
    let y_values = &json.y_values;
    let field = json.field;
    let (coefficients, steps) =
        multilinear_view_helpers::multilinear_interpolate_over_boolean_hypercube(
            y_values,
            field as u128,
        );
    let response = MultilinearOverBooleanHypercubeResponse {
        coefficients,
        steps,
    };
    HttpResponse::Ok().json(response)
}

#[post("/multivariate_interpolation_over_finite_field/")]
async fn multivariate_interpolation_over_finite_field(
    json: web::Json<MultivariateInterpolationRequest>,
) -> impl Responder {
    let num_of_vars = json.num_of_vars;
    let evaluations_points = &json.evaluation_points;
    let y_values = &json.y_values;
    let field = json.field;
    let terms = multivariate_interpolate_over_finite_field(
        num_of_vars,
        evaluations_points,
        y_values,
        field,
    );
    if terms.is_ok() {
        let response = MultivariateInterpolationResponse {
            terms: terms.unwrap(),
        };
        HttpResponse::Ok().json(response)
    } else {
        let response = BadRequestResponse {
            detail: String::from("invalid inputs"),
        };
        HttpResponse::BadRequest().json(response)
    }
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let port = 8080;
    println!("Server is running");

    HttpServer::new(move || {
        App::new()
            .wrap(
                Cors::default()
                    .allowed_origin("https://www.modcc.xyz")
                    .allowed_origin("http://localhost:3000")
                    .allowed_methods(vec!["GET", "POST"])
                    .allowed_headers(vec![http::header::AUTHORIZATION, http::header::ACCEPT])
                    .allowed_header(http::header::CONTENT_TYPE)
                    .max_age(3600),
            )
            .service(lagrange_interpolation_over_ff)
            .service(multilinear_interpolation_over_boolean_hypercube)
            .service(multivariate_interpolation_over_finite_field)
            .service(evaluate_univariate_poly)
    })
    .bind(("127.0.0.1", port))?
    .workers(2)
    .run()
    .await
}
