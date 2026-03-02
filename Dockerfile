# ── Stage 1: Build ────────────────────────────────────────────────────────────
FROM rust:1.85-bookworm AS builder

WORKDIR /build

# Cache dependency build: copy manifests first, build a dummy, then overlay src
COPY Cargo.toml Cargo.lock ./
COPY crates/needletail-core/Cargo.toml crates/needletail-core/Cargo.toml
COPY crates/needletail-py/Cargo.toml   crates/needletail-py/Cargo.toml
COPY crates/needletail-server/Cargo.toml crates/needletail-server/Cargo.toml

# Create dummy source files so cargo can resolve the workspace
RUN mkdir -p crates/needletail-core/src && echo "" > crates/needletail-core/src/lib.rs \
 && mkdir -p crates/needletail-py/src   && echo "" > crates/needletail-py/src/lib.rs \
 && mkdir -p crates/needletail-server/src && echo "fn main() {}" > crates/needletail-server/src/main.rs

# Pre-build dependencies (cached unless Cargo.toml/lock change)
RUN cargo build --release -p needletail-server 2>/dev/null || true

# Copy real source and embedded assets
COPY crates/ crates/
COPY presets/ presets/

# Touch source files to invalidate the dummy build cache
RUN touch crates/needletail-core/src/lib.rs \
          crates/needletail-server/src/main.rs

# Full release build
RUN cargo build --release -p needletail-server

# ── Stage 2: Runtime ─────────────────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/needletail-server /usr/local/bin/

ENV BIND_ADDR=0.0.0.0:8002

EXPOSE 8002

CMD ["needletail-server"]
