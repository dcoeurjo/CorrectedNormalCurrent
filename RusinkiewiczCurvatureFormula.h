/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file RusinkiewiczCurvatureFormula.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2020/01/01
 *
 * Header file for module RusinkiewiczCurvatureFormula.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(RusinkiewiczCurvatureFormula_RECURSES)
#error Recursive header files inclusion detected in RusinkiewiczCurvatureFormula.h
#else // defined(RusinkiewiczCurvatureFormula_RECURSES)
/** Prevents recursive inclusion of headers. */
#define RusinkiewiczCurvatureFormula_RECURSES

#if !defined RusinkiewiczCurvatureFormula_h
/** Prevents repeated inclusion of headers. */
#define RusinkiewiczCurvatureFormula_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <map>
#include "DGtal/base/Common.h"
#include "DGtal/math/linalg/SimpleMatrix.h"
#include "DGtal/math/linalg/EigenDecomposition.h"

#include "geometrycentral/utilities/vector3.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class RusinkiewiczCurvatureFormula
  /**
     Description of class 'RusinkiewiczCurvatureFormula' <p> \brief
     Aim: A helper class that provides static methods to compute
     Rusinkiewicz formulas of curvatures.

     @tparam TRealPoint any model of 3D RealPoint.
     @tparam TRealVector any model of 3D RealVector.
  */
  struct RusinkiewiczCurvatureFormula
  {
    typedef geometrycentral::Vector3            RealPoint;
    typedef geometrycentral::Vector3            RealVector;
    typedef double                         Scalar;
    typedef std::vector< RealPoint >       RealPoints;
    typedef std::vector< RealVector >      RealVectors;
    typedef SimpleMatrix< Scalar, 3, 3 >   RealTensor;
    typedef SimpleMatrix< Scalar, 2, 2 >   RealTensor2D;
    typedef typename RealTensor2D::ColumnVector ColumnVector2D;
    typedef typename RealTensor2D::RowVector    RowVector2D;
    typedef std::size_t                    Index;
    static const Dimension dimension = 3;
    typedef SimpleMatrix< Scalar, 6, 3 >   LSMatrix;
    typedef typename LSMatrix::ColumnVector LSColumnVector;

    //-------------------------------------------------------------------------
  public:
    /// @name Formulas for curvature
    /// @{

    /// Computes a local orthonormal basis on triangle abc.
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @return a pair (u,v) of unit orthogonal vector lying in the plane abc.
    static
    std::pair< RealVector, RealVector > basis
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c )
    {
      RealVector u = ( b - a ).normalize();
      RealVector v = cross(normal( a, b, c ),  u );
      return std::make_pair( u, v );
    }    

    /// Computes a local orthonormal basis on face \a x.
    /// @param x the (ccw ordered) points forming the vertices of a polygonal face.
    /// @return a pair (u,v) of unit orthogonal vector lying in the plane of the face.
    static
    std::pair< RealVector, RealVector > basis
    ( const RealPoints& x )
    {
      if ( x.size() == 3 ) return basis( x[ 0 ], x[ 1 ], x[ 2 ] );
      else                 return basis( barycenter( x ), x[ 0 ], x[ 1 ] );
    }
    
    
    /// Computes the local second fundamental form on triangle abc
    /// given normal vectors \a ua, \a \ub, \a uc at its vertices.
    ///
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @param ua the corrected normal vector at point a
    /// @param ub the corrected normal vector at point b
    /// @param uc the corrected normal vector at point c
    /// @return the 2x2 tensor II according to Rusinkiewicz's formula.
    static
    RealTensor2D secondFundamentalForm
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c,
      const RealVector& ua, const RealVector& ub, const RealVector& uc )
    {
      LSMatrix M;
      LSColumnVector Y;
      prepareFit( M, Y, a, b, c, ua, ub, uc );
      RealTensor2D II;
      bool ok = computeFit( II, M, Y );
      return II;
    }    

    /// Computes the mean curvature on triangle abc
    /// given normal vectors \a ua, \a \ub, \a uc at its vertices.
    ///
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @param ua the corrected normal vector at point a
    /// @param ub the corrected normal vector at point b
    /// @param uc the corrected normal vector at point c
    /// @return the mean curvature according to Rusinkiewicz's formula.
    static
    Scalar meanCurvature
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c,
      const RealVector& ua, const RealVector& ub, const RealVector& uc )
    {
      const auto II = secondFundamentalForm( a, b, c, ua, ub, uc );
      return 0.5 * ( II(0,0)+II(1,1) );
    }

    /// Computes the Gaussian curvature on triangle abc
    /// given normal vectors \a ua, \a \ub, \a uc at its vertices.
    ///
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @param ua the corrected normal vector at point a
    /// @param ub the corrected normal vector at point b
    /// @param uc the corrected normal vector at point c
    /// @return the Gaussian curvature according to Rusinkiewicz's formula.
    static
    Scalar gaussianCurvature
    ( const RealPoint& a, const RealPoint& b, const RealPoint& c,
      const RealVector& ua, const RealVector& ub, const RealVector& uc )
    {
      const auto II = secondFundamentalForm( a, b, c, ua, ub, uc );
      return II(0,0) * II(1,1) - II(0,1) * II(1,0);
    }

    /// Computes second fundamental form of polygonal face \a x given 
    ///  normal vectors at vertices \a u.
    /// @param x the (ccw ordered) points forming the vertices of a polygonal face.
    /// @param u the (ccw ordered) normal vectors at the corresponding vertices in \a x.
    ///
    /// @return the second fundamental form of polygonal face \a x,
    /// expressed in the basis of the triangle (b,x[0],x[1]), or
    /// (x[0], x[1], x[2]) if it is a triangle.
    static
    RealTensor2D secondFundamentalForm( const RealPoints& x, const RealVectors& u )
    {
      ASSERT( x.size() == u.size() );
      if ( x.size() <  3 ) return RealTensor2D();
      if ( x.size() == 3 )
	return secondFundamentalForm( x[ 0 ], x[ 1 ], x[ 2 ],
				      u[ 0 ], u[ 1 ], u[ 2 ] );
      const RealPoint   b = barycenter( x );
      const RealVector ub = averageUnitVector( u );

      const auto   basisp = basis( b, x[ 0 ], x[ 1 ] );
      const auto       up = basisp.first;
      const auto       vp = basisp.second;
      RealTensor2D     II;
      Scalar           w = 0.0;
      for ( Index i = 0; i < x.size(); i++ )
	{
	  const auto basisf = basis( b, x[ i ], x[ (i+1)%x.size() ] );
	  const auto II_i = secondFundamentalForm( b, x[ i ], x[ (i+1)%x.size() ],
						   ub,u[ i ], u[ (i+1)%x.size() ] );
	  const auto w_i  = area( b, x[ i ], x[ (i+1)%x.size() ] );
	  w  += w_i;
	  II += w_i * transformTensor2D( II_i, basisf.first, basisf.second, up, vp );
	}
      return (1.0 / w) * II;
    }

    /// Computes mean curvature of polygonal face \a pts given normal
    /// vectors at vertices \a u.
    /// @param pts the (ccw ordered) points forming the vertices of a polygonal face.
    /// @param u the (ccw ordered) normal vectors at the corresponding vertices in \a pts.
    /// @return the mean curvature of the given polygonal face.
    static
    Scalar meanCurvature( const RealPoints& pts, const RealVectors& u )
    {
      const auto II = secondFundamentalForm( pts, u );
      return 0.5 * ( II(0,0)+II(1,1) );
    }

    /// Computes Gaussian curvature of polygonal face \a pts given normal
    /// vectors at vertices \a u.
    /// @param pts the (ccw ordered) points forming the vertices of a polygonal face.
    /// @param u the (ccw ordered) normal vectors at the corresponding vertices in \a pts.
    /// @return the Gaussian curvature of the given polygonal face.
    static
    Scalar gaussianCurvature( const RealPoints& pts, const RealVectors& u )
    {
      const auto II = secondFundamentalForm( pts, u );
      return II(0,0) * II(1,1) - II(0,1) * II(1,0);
    }

    /// @}

    //-------------------------------------------------------------------------
  public:
    /// @name Least squares services
    /// @{

    static
    bool computeFit( RealTensor2D& II, const LSMatrix& M, const LSColumnVector& Y )
    {
      const auto  tM = M.transpose(); // M(6x3), tM(3*6)
      const auto tMM = tM * M;        // tMM(3x3)
      if ( tMM.determinant() != 0 ) {
	const auto R = tMM.inverse() * ( tM * Y );
	II.setComponent( 0, 0, R[ 0 ] );
	II.setComponent( 0, 1, R[ 1 ] );
	II.setComponent( 1, 0, R[ 1 ] );
	II.setComponent( 1, 1, R[ 2 ] );
	return true;
      }
      return false;
    }
    
    static
    void prepareFit( LSMatrix& M, LSColumnVector& Y,
		     const RealPoint& a, const RealPoint& b, const RealPoint& c,
		     const RealVector& ua, const RealVector& ub, const RealVector& uc )
    {
      std::array< RealVector, 3 >  e = { c - b, a - c, b - a };
      std::array< RealVector, 3 > dn = { uc - ub, ua - uc, ub - ua };
      const auto bf = basis( a, b, c );
      const auto u  = bf.first;
      const auto v  = bf.second;
      for ( Dimension i = 0; i < 3; ++i )
	{
	  Y[ 2*i   ] = dot(dn[ i ], u );
	  Y[ 2*i+1 ] = dot(dn[ i ], v );
	  M.setComponent( 2*i, 0, dot(e[ i ], u ) );
	  M.setComponent( 2*i, 1, dot(e[ i ], v ) );
	  M.setComponent( 2*i, 2, 0.0 );
	  M.setComponent( 2*i+1, 0, 0.0 );
	  M.setComponent( 2*i+1, 1, dot(e[ i ], u ) );
	  M.setComponent( 2*i+1, 2, dot(e[ i ], v ) );
	}
    }
    
    //-------------------------------------------------------------------------
  public:
    /// @name Other geometric services
    /// @{
    
    /// Given a vector of points, returns its barycenter.
    /// @param pts any vector of points
    /// @return the barycenter of these points.
    static 
    RealPoint barycenter( const RealPoints& pts )
    {
      RealPoint b;
      for ( auto p : pts ) b += p;
      b /= pts.size();
      return b;
    }

    /// Computes a unit normal vector to triangle abc
    ///
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @return the unit normal vector to abc, ( ab x ac ) / || ab x ac ||.
    static
    RealVector normal( const RealPoint& a, const RealPoint& b, const RealPoint& c )
    {
      return ( cross(( b - a ),( c - a ) ).normalize() );
    }    

    /// Computes triangle area
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @return the area of triangle abc
    static
    Scalar area( const RealPoint& a, const RealPoint& b, const RealPoint& c )
    {
      return 0.5 * ( cross(( b - a ),( c - a ) ).norm());
    }    

    /// Given a vector of unit vectors, returns their average unit vector.
    /// @param pts any vector of vectors.
    /// @return the average unit vector.
    static 
    RealVector averageUnitVector( const RealVectors& vecs )
    {
      RealVector avg;
      for ( auto v : vecs ) avg += v;
      auto avg_norm = avg.norm();
      return avg_norm != 0.0 ? avg / avg_norm : avg;
    }
    
    /// Transforms tensor \a T from frame (uf,vf) to frame (up,vp).
    static 
    RealTensor2D transformTensor2D( const RealTensor2D& T,
				    const RealVector& uf, const RealVector& vf, 
				    const RealVector& up, const RealVector& vp )
    {
      RealTensor2D U;
      const auto up_uf = dot(up , uf );
      const auto up_vf = dot(up , vf );
      const auto vp_uf = dot(vp , uf );
      const auto vp_vf = dot(vp , vf );
      U.setComponent
	( 0, 0, RowVector2D{ up_uf, up_vf }.dot(T * ColumnVector2D{ up_uf, up_vf }) );
      U.setComponent
	( 0, 1, RowVector2D{ up_uf, up_vf }.dot(T * ColumnVector2D{ vp_uf, vp_vf }) );
      U.setComponent
	( 1, 0, RowVector2D{ vp_uf, vp_vf }.dot(T * ColumnVector2D{ up_uf, up_vf }) );
      U.setComponent
	( 1, 1, RowVector2D{ vp_uf, vp_vf }.dot(T * ColumnVector2D{ vp_uf, vp_vf }) );
      return U;
    }
    
    /// @}
    
  };

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "RusinkiewiczCurvatureFormula.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined RusinkiewiczCurvatureFormula_h

#undef RusinkiewiczCurvatureFormula_RECURSES
#endif // else defined(RusinkiewiczCurvatureFormula_RECURSES)
