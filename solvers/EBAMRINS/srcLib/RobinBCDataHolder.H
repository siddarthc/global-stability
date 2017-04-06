/*
*
*
*
*
*/

#ifndef _ROBINBCDATAHOLDER_H_
#define _ROBINBCDATAHOLDER_H_

#include "LevelData.H"
#include "EBFluxFAB.H"

struct RobinBCData
{
  const Vector<LevelData<EBFluxFAB>* >* dataPtr;
  std::map<Real, int> levelDxMap;
  Real Re;
  Vector<Real> dxVect;
};

class RobinPoissonBCFunctions
{
public:

  static Real extrapToBoundary(const FaceIndex&      a_bndryFace,
                               const Side::LoHiSide& a_side,
                               const int&            a_idir,
                               const EBGraph&        a_ebGraph,
                               const EBFaceFAB&      a_faceData,
                               const int&            a_comp,
                               const RobinBCData&    a_RobBCData,
                               const Real&           a_dx,
                               const DataIndex&      a_dit)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx)->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx)->second : 0;

    Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](a_bndryFace, 0);
    baseDataVal *= a_RobBCData.Re*a_dx*sign(a_side);

    Real extrapValue = -1.e99;
    Side::LoHiSide flipSide = flip(a_side);
    const VolIndex & closeVoF = a_bndryFace.getVoF(flipSide);
    Vector<FaceIndex> nearFaces = a_ebGraph.getFaces(closeVoF, a_idir, flipSide);
    bool hasNearFace = ((nearFaces.size() == 1)  && !nearFaces[0].isBoundary());
    bool hasFarFace  = false;
    bool hasFarFarFace  = false;
    Real nearVal = 0.0;
    Real farVal  = 0.0;
    Real farFarVal  = 0.0;
    if (hasNearFace)
    {
      nearVal = a_faceData(nearFaces[0], 0);
      const VolIndex & nextVoF = nearFaces[0].getVoF(flipSide);
      Vector<FaceIndex> farFaces = a_ebGraph.getFaces(nextVoF, a_idir, flipSide);
      hasFarFace = ((farFaces.size() == 1) && !farFaces[0].isBoundary());
      if (hasFarFace)
      {
        farVal = a_faceData(farFaces[0], 0);
        const VolIndex & nextNextVoF = farFaces[0].getVoF(flipSide);
        Vector<FaceIndex> farFarFaces = a_ebGraph.getFaces(nextNextVoF, a_idir, flipSide);
        hasFarFarFace = ((farFarFaces.size() == 1) && !farFarFaces[0].isBoundary());
        if (hasFarFarFace)
        {
          farFarVal = a_faceData(farFarFaces[0], 0);
        }
      }
    }

  // if (hasNearFace && hasFarFace && hasFarFarFace)
  //   {
  //     extrapValue = 3.*nearVal - 3.*farVal + farFarVal;
  //   }
  // else if (hasNearFace && hasFarFace)
  //   {
  //     extrapValue = 2.*nearVal - farVal;
  //   }
    if (hasNearFace && hasFarFace)
    {
      //quadratic extrapolation on v with the bc dv/dx + Re*U.n*v = 0 applied
      extrapValue = (4.*nearVal - farVal)/(3. - 2.*baseDataVal);
    }
    else if (hasNearFace)
    {
      extrapValue = nearVal/(1. + baseDataVal);
    }
    else
    {
      extrapValue = 0.0; //for want of a better option.
    }
    return extrapValue;
  }

  /*********/
  static void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                          const BaseFab<Real>&  a_phi,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const Side::LoHiSide& a_side,
                          const DataIndex&      a_dit,
                          const Real&           a_time,
                          const bool&           a_useHomogeneous,
                          const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;

    for (int comp=0; comp<a_phi.nComp(); comp++)
    {
      const Box& box = a_faceFlux.box();

      int iside;

      if (a_side == Side::Lo)
        {
          iside = 1;
        }
      else
        {
          iside = -1;
        }

      Real ihdx;
      ihdx = 2.0 / a_dx[a_idir];
      Real weight = iside*ihdx;
      BoxIterator bit(box);

      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        IntVect ivNeigh = iv;
        ivNeigh[a_idir] += sign(a_side);
        const VolIndex vof      = VolIndex(iv,     0);
        const VolIndex vofNeigh = VolIndex(ivNeigh,0);

        const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);

        IntVect baseIV = iv;
        baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
        IntVect baseIVNeigh = baseIV;
        baseIVNeigh[a_idir] += sign(a_side);
        const VolIndex basevof      = VolIndex(baseIV,     0);
        const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);

        const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

        Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
        baseDataVal *= -1.*a_RobBCData.Re;
        Real factor = baseDataVal/(baseDataVal + weight);

        Real phiVal = a_phi(iv,comp);
        
        // Here's how it works: dphi/dn = wb*phib + Sum(wi*phii) = C*phib
        //     => phib = Sum(wi*phii)/(C-wb)
        //     => dphi/dn = C* Sum(wi*phii)/(C-wb)
        a_faceFlux(iv,comp) = factor*weight*phiVal;
      }
    }
  }

  /*********/
  static void getHigherOrderFaceFlux(BaseFab<Real>&        a_faceFlux,
                                     const BaseFab<Real>&  a_phi,
                                     const RealVect&       a_probLo,
                                     const RealVect&       a_dx,
                                     const int&            a_idir,
                                     const Side::LoHiSide& a_side,
                                     const DataIndex&      a_dit,
                                     const Real&           a_time,
                                     const bool&           a_useHomogeneous,
                                     const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;

    for (int comp=0; comp<a_phi.nComp(); comp++)
    {
      const Box& box = a_faceFlux.box();

      int iside;

      if (a_side == Side::Lo)
        {
          iside = 1;
        }
      else
        {
          iside = -1;
        }
    
      // Here's how it works: dphi/dn = wb*phib + Sum(wi*phii) = C*phib
      //     => phib = Sum(wi*phii)/(C-wb)
      //     => dphi/dn = C* Sum(wi*phii)/(C-wb)

      BoxIterator bit(box);

      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        IntVect ivNeigh = iv;
        ivNeigh[a_idir] += sign(a_side);
        const VolIndex vof      = VolIndex(iv,     0);
        const VolIndex vofNeigh = VolIndex(ivNeigh,0);
        const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);

        IntVect baseIV = iv;
        baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
        IntVect baseIVNeigh = baseIV;
        baseIVNeigh[a_idir] += sign(a_side);
        const VolIndex basevof      = VolIndex(baseIV,     0);
        const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);

        const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

        Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
        baseDataVal *= -1.*a_RobBCData.Re;
        Real factor = baseDataVal/(baseDataVal + iside*8./(3.*a_dx[a_idir]));

        Real phi0 = a_phi(iv,comp);
        iv[a_idir] += iside;
        Real phi1 = a_phi(iv,comp);
        iv[a_idir] -= iside;
        a_faceFlux(iv,comp) = factor*-iside*(phi1 - 9.0*phi0)/(3.0*a_dx[a_idir]);
      }
    }
  }

/*********/
  static void getFaceFlux(Real&                 a_faceFlux,
                          const VolIndex&       a_vof,
                          const int&            a_comp,
                          const EBCellFAB&      a_phi,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const Side::LoHiSide& a_side,
                          const DataIndex&      a_dit,
                          const Real&           a_time,
                          const bool&           a_useHomogeneous,
                          const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;
    int iside = -sign(a_side);
    const Real ihdx = 2.0 / a_dx[a_idir];

    const EBISBox& ebisBox = a_phi.getEBISBox();
    const ProblemDomain& domainBox = ebisBox.getDomain();
    a_faceFlux = 0.0;
    Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    if (faces.size() > 0)
    {
      if (faces.size()==1)
      {
        IntVectSet cfivs;
        FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                         cfivs,
                                                         ebisBox,
                                                         domainBox);
        for (int isten=0; isten < faceSten.size(); isten++)
        {
          const Real& weight = faceSten.weight(isten);
          const FaceIndex& face = faceSten.face(isten);
          Real thisFaceFlux;
          const RealVect centroid = ebisBox.centroid(face);

          IntVect baseIV = face.gridIndex(a_side);
          baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
          IntVect baseIVNeigh = baseIV;
          baseIVNeigh[a_idir] += sign(a_side);
          const VolIndex basevof      = VolIndex(baseIV,     0);
          const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);

          const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

          Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
          baseDataVal *= -1.*a_RobBCData.Re;
          Real factor = baseDataVal/(baseDataVal + iside*ihdx);

          thisFaceFlux = factor*iside*ihdx*a_phi(a_vof,a_comp);

          a_faceFlux += thisFaceFlux*weight;
        }
        a_faceFlux *= ebisBox.areaFrac(faces[0]);
      }
      else
      {
        MayDay::Error("RobinBCFunctions::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
      }
    }
  }

/*********/
  static void getHigherOrderFaceFlux(Real&                 a_faceFlux,
                                     const VolIndex&       a_vof,
                                     const int&            a_comp,
                                     const EBCellFAB&      a_phi,
                                     const RealVect&       a_probLo,
                                     const RealVect&       a_dx,
                                     const int&            a_idir,
                                     const Side::LoHiSide& a_side,
                                     const DataIndex&      a_dit,
                                     const Real&           a_time,
                                     const bool&           a_useHomogeneous,
                                     const RobinBCData&    a_RobBCData)
  {
    const EBISBox& ebisBox = a_phi.getEBISBox();
    const ProblemDomain& domainBox = ebisBox.getDomain();
    a_faceFlux = 0.0;
    Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    if (faces.size() > 0)
    {
      if (faces.size()==1)
      {
        IntVectSet cfivs;
        FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                         cfivs,
                                                         ebisBox,
                                                         domainBox);
        for (int isten=0; isten < faceSten.size(); isten++)
        {
          const Real& weight = faceSten.weight(isten);
          const FaceIndex& face = faceSten.face(isten);
          Real thisFaceFlux;
          const RealVect centroid = ebisBox.centroid(face);
          getHigherOrderFaceFlux(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                 a_side,a_dit,a_time,false,centroid,a_useHomogeneous,a_RobBCData);

          a_faceFlux += thisFaceFlux*weight;
        }
        a_faceFlux *= ebisBox.areaFrac(faces[0]);
      }
      else
      {
        MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
      }
    }
  } 

/*********/
  static void getHigherOrderFaceFlux(Real&                 a_faceFlux,
                                     const FaceIndex&      a_face,
                                     const int&            a_comp,
                                     const EBCellFAB&      a_phi,
                                     const RealVect&       a_probLo,
                                     const RealVect&       a_dx,
                                     const int&            a_idir,
                                     const Side::LoHiSide& a_side,
                                     const DataIndex&      a_dit,
                                     const Real&           a_time,
                                     const bool&           a_useAreaFrac,
                                     const RealVect&       a_centroid,
                                     const bool&           a_useHomogeneous,
                                     const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;
    int iside = -sign(a_side);
    const Real ihdx = 2.0 / a_dx[a_idir];

    const EBISBox& ebisBox = a_phi.getEBISBox();
    const VolIndex& vof = a_face.getVoF(flip(a_side));

    IntVect baseIV = a_face.gridIndex(flip(a_side));
    baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
    IntVect baseIVNeigh = baseIV;
    baseIVNeigh[a_idir] += sign(a_side);
    const VolIndex basevof      = VolIndex(baseIV,     0);
    const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);
    const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

    Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
    baseDataVal *= -1.*a_RobBCData.Re;

    //higher-order Dirichlet bc
    Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
    if (facesInsideDomain.size() == 1)
    {
      Real factor = baseDataVal/(baseDataVal + iside*8./(3.*a_dx[a_idir]));
 
      const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
      a_faceFlux = factor * iside * (9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp))/(3.0*a_dx[a_idir]);
    }
    else
    {
      Real factor = baseDataVal/(baseDataVal + iside*ihdx);
      a_faceFlux = factor * iside * ihdx * a_phi(vof,a_comp);
    }
  }

/*********/
  static void getInhomFaceFlux(Real&                 a_faceFlux,
                               const VolIndex&       a_vof,
                               const int&            a_comp,
                               const EBCellFAB&      a_phi,
                               const RealVect&       a_probLo,
                               const RealVect&       a_dx,
                               const int&            a_idir,
                               const Side::LoHiSide& a_side,
                               const DataIndex&      a_dit,
                               const Real&           a_time,
                               const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);
  
    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;
    int iside = -sign(a_side);
    const Real ihdx = 2.0 / a_dx[a_idir];
  
    const EBISBox& ebisBox = a_phi.getEBISBox();
    const ProblemDomain& domainBox = ebisBox.getDomain();
    a_faceFlux = 0.0;
    Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    if (faces.size() > 0)
    {
      if (faces.size()==1)
      {
        IntVectSet cfivs;
        FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                         cfivs,
                                                         ebisBox,
                                                         domainBox);
        for (int isten=0; isten < faceSten.size(); isten++)
        {
          const Real& weight = faceSten.weight(isten);
          const FaceIndex& face = faceSten.face(isten);
          Real thisFaceFlux;
          const RealVect centroid = ebisBox.centroid(face);

          IntVect baseIV = face.gridIndex(a_side);
          baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
          IntVect baseIVNeigh = baseIV;
          baseIVNeigh[a_idir] += sign(a_side);
          const VolIndex basevof      = VolIndex(baseIV,     0);
          const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);

          const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

          Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
          baseDataVal *= -1.*a_RobBCData.Re;
          Real factor = baseDataVal/(baseDataVal + iside*ihdx);
          thisFaceFlux = factor*iside*ihdx*a_phi(a_vof,a_comp);
          a_faceFlux += thisFaceFlux*weight;
        }
        a_faceFlux *= ebisBox.areaFrac(faces[0]);
      }
      else
      {
        MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
      }
    }
  }

/*********/
  static void getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
                                          const VolIndex&       a_vof,
                                          const int&            a_comp,
                                          const EBCellFAB&      a_phi,
                                          const RealVect&       a_probLo,
                                          const RealVect&       a_dx,
                                          const int&            a_idir,
                                          const Side::LoHiSide& a_side,
                                          const DataIndex&      a_dit,
                                          const Real&           a_time,
                                          const RobinBCData&    a_RobBCData)
  {
    const EBISBox& ebisBox = a_phi.getEBISBox();
    const ProblemDomain& domainBox = ebisBox.getDomain();
    a_faceFlux = 0.0;
    Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
    if (faces.size() > 0)
    {
      if (faces.size()==1)
      {
        IntVectSet cfivs;
        FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                         cfivs,
                                                         ebisBox,
                                                         domainBox);
        for (int isten=0; isten < faceSten.size(); isten++)
        {
          const Real& weight = faceSten.weight(isten);
          const FaceIndex& face = faceSten.face(isten);
          Real thisFaceFlux;
          const RealVect centroid = ebisBox.centroid(face);
          getHigherOrderInhomFaceFlux(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                      a_side,a_dit,a_time,false,centroid,a_RobBCData);

          a_faceFlux += thisFaceFlux*weight;
          }
        a_faceFlux *= ebisBox.areaFrac(faces[0]);
      }
      else
      {
        MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
      }
    }
  }

/*********/
  static void getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
                                          const FaceIndex&      a_face,
                                          const int&            a_comp,
                                          const EBCellFAB&      a_phi,
                                          const RealVect&       a_probLo,
                                          const RealVect&       a_dx,
                                          const int&            a_idir,
                                          const Side::LoHiSide& a_side,
                                          const DataIndex&      a_dit,
                                          const Real&           a_time,
                                          const bool&           a_useAreaFrac,
                                          const RealVect&       a_centroid,
                                          const RobinBCData&    a_RobBCData)
  {
//    int lev = a_RobBCData.levelDxMap.find(a_dx[a_idir])->second;
//    int lev = a_RobBCData.levelDxMap.at(a_dx[a_idir]);

    std::map<Real, int>::const_iterator it = a_RobBCData.levelDxMap.find(a_dx[a_idir]);
    int lev = (it != a_RobBCData.levelDxMap.end()) ? a_RobBCData.levelDxMap.find(a_dx[a_idir])->second : 0;
    int iside = -sign(a_side);
    const Real ihdx = 2.0 / a_dx[a_idir];

    const EBISBox& ebisBox = a_phi.getEBISBox();
    const VolIndex& vof = a_face.getVoF(flip(a_side));

    IntVect baseIV = a_face.gridIndex(flip(a_side));
    baseIV *= a_dx[a_idir]/a_RobBCData.dxVect[lev];
    IntVect baseIVNeigh = baseIV;
    baseIVNeigh[a_idir] += sign(a_side);
    const VolIndex basevof      = VolIndex(baseIV,     0);
    const VolIndex basevofNeigh = VolIndex(baseIVNeigh,0);

    const FaceIndex baseface = FaceIndex(basevof,basevofNeigh,a_idir);

    Real baseDataVal = (*(*a_RobBCData.dataPtr)[lev])[a_dit][a_idir](baseface, a_idir);
    baseDataVal *= -1.*a_RobBCData.Re;

    //higher-order Dirichlet bc
    Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
    if (facesInsideDomain.size() == 1)
    {
      Real factor = baseDataVal/(baseDataVal + iside*8./(3.*a_dx[a_idir]));

      const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
      a_faceFlux = factor * iside * (9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp))/(3.0*a_dx[a_idir]);
    }
    else
    {
      Real factor = baseDataVal/(baseDataVal + iside*ihdx);
      a_faceFlux = factor * iside * ihdx * a_phi(vof,a_comp);
    }
  }
};

#endif
