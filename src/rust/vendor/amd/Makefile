all:
	cargo run --bin simple
	cargo run --bin order
	cargo run --bin jumbled

check:
	cargo clippy -- \
-A clippy::needless_range_loop \
-A clippy::many_single_char_names \
-A clippy::too_many_arguments \
-A clippy::int_plus_one \
-A clippy::comparison_chain

# rustup target add x86_64-unknown-linux-musl
dist:
	cargo build --release --target=x86_64-unknown-linux-musl
