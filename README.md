# Table of Contents

# Introduction

This repository is meant to allow for active subspace sensitivity analysis of plasma dynamics using a particle-in-cell (PIC) method.

# General Usage

The [PIC subdirectory](../+PIC) contains two main functions - one is designed to create movies and images of phase-space dynamics, and the other is for creating the functions used in the sensitivity analysis.  The folders [`Initilization`](../+PIC/Initilization) and [`QOI_calc`](../+PIC/QOI_calc) contain functions for initializing the particle distributions and calculating the quantities of interest, respectively.

The [sensitivity subdirectory](../+Sensitivity) contains functions for gradient-based and linear-fit sensitivity metrics.  Associated functions designed to create plots based on the data are also included.

## Examples

Usage examples for two Landau Damping and one Two-Stream Instability (split into two parts - one with a tight beam and one without) are included.

# Particle-In-Cell Method

Using typical parameter values of ..., each evaluation of the PIC algorithm takes approximately five minutes.  Time complexity of the algorithm is roughly quadratic in number of particles, linear in number of spatial gridpoints, and linear in both temporal gridpoints and maximal value.  It's recommended that when using the [`movie_run.m`](../+PIC/movie_run.m) more accurate particles and spatial gridpoints are used to ensure proper resolution.  Maximal time is also taken in all examples to be close to the minimum possible, and longer time tends to show more clearly the dynamics.

## Functions

[`movie_run( DT, NT, NG, N, distribution, params, saveFrameNum, movieFrameNum)`](../PIC/movie_run.m)
| Parameters | Meaning |
| --- | --- |
| `DT` | Temporal step size |
| `NT` | Number of temporal steps |
| `NG` | Number of spatial gridpoints |
| `N` | Number of computational particles |
| `distribution` | Name of the distribution - for use with functions in [`Initilization`](../+PIC/Initilization) and [`QOI_calc`](../+PIC/QOI_calc) |
| `params` | Parameters for the initilization of the distribution - the spatial length is assumed to come first |
| `saveFrameNum` | Saves a .png file of the phase-space dynamics at every `saveFrameNum` temporal step (must be a multiple of `movieFrameNum`) |
| `movieFrameNum` | Adds a frame to the .mp4 movie of phase-space dynamics at every `movieFrameNum` temporal step |

## Examples

# Active Subspace Sensitivity Methods

## Functions

![Alt text](/relative/path/to/img.jpg?raw=true "title")
<!--- can use a branch with the examples to not clutter the main directory --->
<!--- also potentially try embedding movies? --->

## Examples