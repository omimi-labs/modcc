use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use actix_web::{web, App, HttpResponse, HttpServer, Responder};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct User {
    name: String,
}

type UserDb = Arc<Mutex<HashMap<u32, User>>>;

#[actix_web::get("/greet/{id}")]
async fn greet(user_id: web::Path<u32>) -> impl Responder {
    format!("Hello World {}", user_id)
}

#[actix_web::post("/user/")]
async fn create_user(user_data: web::Json<User>, db: web::Data<UserDb>) -> impl Responder {
    let mut db = db.lock().unwrap();
    let new_id = db.keys().max().unwrap_or(&0) + 1;
    let name = user_data.name.clone();
    db.insert(new_id, user_data.into_inner());
    HttpResponse::Created().json(User { name })
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let port = 8080;
    println!("Server is running");

    let user_db: UserDb = Arc::new(Mutex::new(HashMap::<u32, User>::new()));

    HttpServer::new(move || {
        let app_data = web::Data::new(user_db.clone());
        App::new()
            .app_data(app_data)
            .service(greet)
            .service(create_user)
    })
    .bind(("127.0.0.1", port))?
    .workers(2)
    .run()
    .await
}
