# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [5.2.0]
### Changed
- Require Julia 1.10 or later.
- Optimize nucleotide complementation.

## [5.1.4]
### Fixed
- Qualify the `Char` constructor extension as `Base.Char`.

## [5.1.3]
### Changed
- Replace SnoopPrecompile with PrecompileTools.

### Fixed
- Improve errors when converting an invalid character to a predefined biological symbol.

## [5.1.2]
### Added
- Define `gap(Char)` here so downstream packages do not need to extend it.
- Add precompilation with SnoopPrecompile.

## [5.1.1]
### Changed
- Optimize symbol hashing.
- Add a workflow for testing downstream packages.

## [5.1.0]
### Added
- Add `stringbyte` for symbols representable as ASCII characters.

## [5.0.0]
### Added
- Make biological symbols broadcast as scalars.

### Changed
- Store small lookup tables as tuples and improve type inference.

### Removed
- Remove iteration, integer conversion, and arithmetic on biological symbols.
- Remove Automa-based parsing.

## [4.0.4]
### Changed
- Update Julia compatibility for long-term support releases.

## [4.0.3]
- Update dependencies

## [4.0.2]
- Changed lookup tables to Tuple for performance improvements

## [4.0.1]
- Fixed missing iterate and length methods

## [4.0.0]
### Added
- An abstract type `BioSymbol` has been added.

### Changed
- `tryparse` now returns `nothing` instead of a `Nullable`.
- The project now uses Project TOML files instead of a REQUIRE file.
- Documentation has been changed to use Documenter.jl native html generation.
- Documentation has a new layout with a manual section, and a library section.

## [3.0.2]
### Added
- Conversion between the DNA and RNA symbol types.

## [3.0.1]
### Added
- Some additional compatiblity changes for julia v0.7 / v1.0

## [3.0.0]
### Added
- Support for julia v0.7 / v1.0. 

### Removed
- :exclamation: Dropped support for julia v0.6.

## [3.0.0]
### Added
- Support for julia v0.7 / v1.0. 

### Removed
- :exclamation: Dropped support for julia v0.6.

## [2.0.1] - 2018-07-02
### Added
- Read and write methods to fix errors when using `pmap`.
- Contributing files and GitHub templates. 

## [2.0.0] - 2018-03-13
### Added
- Support for julia v0.7.

### Removed
- :exclamation: Dropped support for julia v0.5.

## [1.2.0] - 2017-07-26
### Added
- Automa based amino acid parser.

### Changed
- Documentation structure.
- Test coverage increases.

## [1.1.1] - 2017-07-04
### Changed
- Documentation changes.
- Generalised the `show` method for `NucleicAcid` types.


[Unreleased]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.2.0...HEAD
[5.2.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.1.4...v5.2.0
[5.1.4]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.1.3...v5.1.4
[5.1.3]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.1.2...v5.1.3
[5.1.2]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.1.1...v5.1.2
[5.1.1]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.1.0...v5.1.1
[5.1.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v5.0.0...v5.1.0
[5.0.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v4.0.4...v5.0.0
[4.0.4]: https://github.com/BioJulia/BioSymbols.jl/compare/v4.0.3...v4.0.4
[4.0.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v3.0.2...v4.0.0
[3.0.2]: https://github.com/BioJulia/BioSymbols.jl/compare/v3.0.1...v3.0.2
[3.0.1]: https://github.com/BioJulia/BioSymbols.jl/compare/v3.0.0...v3.0.1
[3.0.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v2.0.1...v3.0.0
[2.0.1]: https://github.com/BioJulia/BioSymbols.jl/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v1.2.0...v2.0.0
[1.2.0]: https://github.com/BioJulia/BioSymbols.jl/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/BioJulia/BioSymbols.jl/compare/v1.1.0...v1.1.1
