/**
 Copyright (c) 2020,  CNRS
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,
 David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 LIRIS (CNRS, UMR 5205), CNRS, France
 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the <organization> nor the names of its 
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT
HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once
#include <iostream>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <array>
#include <Eigen/Dense>
#include <assert.h>
/**
 * @file CorrectedNormalCurrentEigen.h
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205), CNRS, France
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2020/12/09
 *
 */

/**
 * This header file contains some stand-alone CorrectedNormalCurrent formulas for triangles,
 * using eigen as linear algebra backend.
 *
 * For additional formulas (without the normal vector interpolations,
 * for convex polygonal faces...) and detailed documentation,
 * please refer to CorrectedNormalCurrent.h and SphericalTriangle.h
 *
 */

namespace CorrectedNormalCurrentEigen {
  struct SphericalTriangle {
    
    ///Spherical point data type
    typedef Eigen::Vector3d Vector3;
    
    static
    bool isDegenerate(const Vector3 &a, const Vector3 &b, const Vector3 &c)
    {
      double d[ 3 ] = { ( a - b ).norm(),( a - c ).norm(),( b - c ).norm() };
      // Checks that the spherical triangle is small or thin.
      if ( ( d[ 0 ] < 1e-8 ) || ( d[ 1 ] < 1e-8 ) || ( d[ 2 ] < 1e-8 ) )
        return true;
      // Checks that the spherical triangle is flat.
      size_t m = 0;
      if ( d[ 1 ] > d[ m ] ) m = 1;
      if ( d[ 2 ] > d[ m ] ) m = 2;
      return ( fabs( d[ m ] - d[ (m+1)%3 ] - d[ (m+2)%3 ] ) < 1e-8 );
    }
    
    /// @return the polar triangle associated with this triangle.
    static
    void polarTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c,
                       Vector3& Ap, Vector3& Bp, Vector3& Cp)
    {
      Ap = b.cross(c);
      Bp = c.cross(a);
      Cp = a.cross(b);
      // Reorient points.
      if ( Ap.dot( a ) < 0.0 ) Ap = -Ap;
      if ( Bp.dot( b ) < 0.0 ) Bp = -Bp;
      if ( Cp.dot( c ) < 0.0 ) Cp = -Cp;
    }
    
    /// Returns the interior angles of the spherical triangle ABC.
    /// @param[out] alpha the interior angle at vertex A.
    /// @param[out] beta  the interior angle at vertex B.
    /// @param[out] gamma the interior angle at vertex C.
    static
    void interiorAngles( const Vector3 &a, const Vector3 &b, const Vector3 &c,
                        double& alpha, double& beta, double& gamma )
    {
      Vector3 Ta,Tb,Tc;
      polarTriangle(a,b,c,Ta,Tb,Tc);
      Ta /= Ta.norm();
      Tb /= Tb.norm();
      Tc /= Tc.norm();
      if ( Ta == Vector3::Zero() || Tb == Vector3::Zero() || Tc == Vector3::Zero() )
        alpha = beta = gamma = 0.0;
      else
      {
        double ca = std::max( -1.0, std::min( 1.0, Tb.dot( Tc ) ) );
        double cb = std::max( -1.0, std::min( 1.0, Tc.dot( Ta ) ) );
        double cc = std::max( -1.0, std::min( 1.0, Ta.dot( Tb ) ) );
        alpha     = acos( ca );
        beta      = acos( cb );
        gamma     = acos( cc );
      }
    }
    
    /// @return the (unsigned) area of the spherical triangle (below 2pi).
    static
    double area(const Vector3 &a, const Vector3 &b, const Vector3 &c)
    {
      double alpha, beta, gamma;
      if ( isDegenerate(a,b,c) ) return 0.0;
      interiorAngles( a, b, c, alpha, beta, gamma );
      return ( (alpha == 0.0) || (beta == 0.0) || (gamma == 0.0) )
      ? 0.0 : 2.0*M_PI - alpha - beta - gamma;
    }
    
    /// @return the (signed) area of the spherical triangle (below 2pi).
    static
    double algebraicArea(const Vector3 &a, const Vector3 &b, const Vector3 &c) 
    {
      double  S = area(a,b,c);
      Vector3 M = a + b + c;
      Vector3 X = ( b - a ).cross( c - a );
      if ( M.lpNorm<1>() <= 1e-8 || X.lpNorm<1>() <= 1e-8 ) return 0.0;
      return M.dot( X ) < 0.0 ? -S : S;
    }
  };
  
  
  ///---------------------- Main functions ----------------
  
  /// Computes mu0 measure (area) of triangle abc given an interpolated
  /// corrected normal vector \a ua, \a \ub, \a uc.
  /// @param a any point
  /// @param b any point
  /// @param c any point
  /// @param ua the corrected normal vector at point a
  /// @param ub the corrected normal vector at point b
  /// @param uc the corrected normal vector at point c
  /// @param unit_u when 'true' considers that interpolated
  /// corrected normals should be made unitary, otherwise
  /// interpolated corrected normals may have smaller norms.
  /// @return the mu0-measure of triangle abc, i.e. its area.
  static
  double mu0InterpolatedU( const Eigen::Vector3d & a, const Eigen::Vector3d & b, const Eigen::Vector3d & c,
                         const Eigen::Vector3d & ua, const Eigen::Vector3d & ub, const Eigen::Vector3d & uc,
                         bool unit_u = false)
  {
    // MU0=1/2*det( uM, B-A, C-A )
    //    =  1/2 < ( (u_A + u_B + u_C)/3.0 ) | (AB x AC ) >
    Eigen::Vector3d uM = ( ua+ub+uc ) / 3.0;
    if ( unit_u )
    {
      auto uM_norm = uM.norm();
      uM = uM_norm == 0.0 ? uM : uM / uM_norm;
    }
    return 0.5 * (( b - a ).cross( c - a )).dot( uM );
  }
  
  
  /// Computes mu1 measure (mean curvature) of triangle abc given an interpolated
  /// corrected normal vector \a ua, \a \ub, \a uc.
  /// @param a any point
  /// @param b any point
  /// @param c any point
  /// @param ua the corrected normal vector at point a
  /// @param ub the corrected normal vector at point b
  /// @param uc the corrected normal vector at point c
  /// @param unit_u when 'true' considers that interpolated
  /// corrected normals should be made unitary, otherwise
  /// interpolated corrected normals may have smaller norms.
  /// @return the mu1-measure of triangle abc, i.e. its mean curvature.
  double mu1InterpolatedU( const Eigen::Vector3d& a, const Eigen::Vector3d & b, const Eigen::Vector3d& c,
                          const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
                          bool unit_u = false)
  {
    // MU1=1/2( | uM u_C-u_B A | + | uM u_A-u_C B | + | uM u_B-u_A C |
    Eigen::Vector3d uM = ( ua+ub+uc ) / 3.0;
    if ( unit_u ) uM /= uM.norm();
    return 0.25 * ( uM.cross( uc - ub ).dot( a )
                   + uM.cross( ua - uc ).dot( b )
                   + uM.cross( ub - ua ).dot( c ) );
  }
  
  
  /// Computes mu2 measure (Gaussian curvature) of triangle abc given an interpolated
  /// corrected normal vector \a ua, \a \ub, \a uc.
  /// @param a any point
  /// @param b any point
  /// @param c any point
  /// @param ua the corrected normal vector at point a
  /// @param ub the corrected normal vector at point b
  /// @param uc the corrected normal vector at point c
  /// @param unit_u when 'true' considers that interpolated
  /// corrected normals should be made unitary, otherwise
  /// interpolated corrected normals may have smaller norms.
  /// @return the mu2-measure of triangle abc, i.e. its Gaussian curvature.
  static
  double mu2InterpolatedU( const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                          const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
                          bool unit_u = false)
  {

    // Using non unitary interpolated normals give
    // MU2=1/2*det( uA, uB, uC )
    // When normals are unitary, it is the area of a spherical triangle.
    if ( unit_u )
      return SphericalTriangle::algebraicArea( ua,ub,uc);
    else
      return 0.5 * ( ua.cross( ub ).dot( uc ) );
  }
  
  
  /// Computes muXY measure (anisotropic curvature) of triangle abc given an interpolated
  /// corrected normal vector \a ua, \a \ub, \a uc.
  /// @param a any point
  /// @param b any point
  /// @param c any point
  /// @param ua the corrected normal vector at point a
  /// @param ub the corrected normal vector at point b
  /// @param uc the corrected normal vector at point c
  /// @return the muXY-measure of triangle abc, i.e. its anisotropic curvature.
  static
  Eigen::Matrix3d muXYInterpolatedU( const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                                   const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
                                    bool unit_u = false)
  {
    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    Eigen::Vector3d uM = ( ua+ub+uc ) / 3.0;
    if ( unit_u ) uM /= uM.norm();
    const Eigen::Vector3d uac = uc - ua;
    const Eigen::Vector3d uab = ub - ua;
    const Eigen::Vector3d  ab = b - a;
    const Eigen::Vector3d  ac = c - a;
    for ( size_t i = 0; i < 3; ++i )
    {
      Eigen::Vector3d X = Eigen::Vector3d::Zero();
      X(i) = 1.0 ;
      for ( size_t j = 0; j < 3; ++j )
      {
        // Since RealVector Y = RealVector::base( j, 1.0 );
        // < Y | uac > = uac[ j ]
        const double tij =
        0.5 * uM.dot( uac[ j ] * X.cross( ab )
                     - uab[ j ] * X.cross( ac ) );
        T(i,j) = tij;
      }
    }
    return T;
  }
  
  
  
  ///---------------------- Helper functions ----------------
  ///
  ///
  /// Computing principal curvatures k1 and k2 from tensor
  /// @param tensor The muXY integrated tensor
  /// @param area Area of the face
  /// @param N the normal vector
  /// @return a pair of principal directions.
  static
  std::pair<Eigen::Vector3d, Eigen::Vector3d>
  curvDirFromTensor(const Eigen::Matrix3d &tensor,
                    const double area,
                    const Eigen::Vector3d &N)
  {
    auto  Mt=tensor.transpose();
    auto M=tensor;
    M += Mt;
    M *= 0.5;
    const double   coef_N = 1000.0 * area;
    // Adding 1000 area n x n to anisotropic measure
    for ( int j = 0; j < 3; j++ )
    for ( int k = 0; k < 3; k++ )
    M( j, k ) += coef_N * N[ j ] * N[ k ];
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
    if (eigensolver.info() != Eigen::Success) abort();
    
    //SelfAdjointEigenSolver returns sorted eigenvalues, no
    //need to reorder the eigenvectors.
    assert(eigensolver.eigenvalues()(0) <= eigensolver.eigenvalues()(1) );
    assert(eigensolver.eigenvalues()(1) <= eigensolver.eigenvalues()(2) );    
    Eigen::Vector3d v1 = eigensolver.eigenvectors().col(1);
    Eigen::Vector3d v2 = eigensolver.eigenvectors().col(0);
    return std::pair<Eigen::Vector3d,Eigen::Vector3d>(v1,v2);
  }
}
