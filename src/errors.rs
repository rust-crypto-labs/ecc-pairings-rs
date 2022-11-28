#[derive(Debug)]
pub enum ErrorKind {
    InvalidInput(&'static str),
    NonQuadraticResidue,
    NoInverse,
}
