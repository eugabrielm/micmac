/*Header-MicMac-eLiSe-25/06/2007

    MicMac : Multi Image Correspondances par Methodes Automatiques de Correlation
    eLiSe  : ELements of an Image Software Environnement

    www.micmac.ign.fr

   
    Copyright : Institut Geographique National
    Author : Marc Pierrot Deseilligny
    Contributors : Gregoire Maillet, Didier Boldo.

[1] M. Pierrot-Deseilligny, N. Paparoditis.
    "A multiresolution and optimization-based image matching approach:
    An application to surface reconstruction from SPOT5-HRS stereo imagery."
    In IAPRS vol XXXVI-1/W41 in ISPRS Workshop On Topographic Mapping From Space
    (With Special Emphasis on Small Satellites), Ankara, Turquie, 02-2006.

[2] M. Pierrot-Deseilligny, "MicMac, un lociel de mise en correspondance
    d'images, adapte au contexte geograhique" to appears in 
    Bulletin d'information de l'Institut Geographique National, 2007.

Francais :

   MicMac est un logiciel de mise en correspondance d'image adapte 
   au contexte de recherche en information geographique. Il s'appuie sur
   la bibliotheque de manipulation d'image eLiSe. Il est distibue sous la
   licences Cecill-B.  Voir en bas de fichier et  http://www.cecill.info.


English :

    MicMac is an open source software specialized in image matching
    for research in geographic information. MicMac is built on the
    eLiSe image library. MicMac is governed by the  "Cecill-B licence".
    See below and http://www.cecill.info.

Header-MicMac-eLiSe-25/06/2007*/
#include "StdAfx.h"

/*************************************************/
/*                                               */
/*          cCC_Appli                            */
/*                                               */
/*************************************************/

cFileOriMnt cAppli_Ortho::GetOriMnt(const std::string & aName) const
{
    return StdGetObjFromFile<cFileOriMnt>
           (
               mWorkDir+aName,
               StdGetFileXMLSpec("ParamChantierPhotogram.xml"),
               "FileOriMnt",
               "FileOriMnt"
           );
}

cFileOriMnt * cAppli_Ortho::GetMtdMNT()
{
   if (mCO.FileMNT().IsInit())
     return new cFileOriMnt(GetOriMnt (mCO.FileMNT().Val()));

   return 0;
}

std::vector<cMetaDataPartiesCachees>  * cAppli_Ortho::VMDPC()
{
   std::vector<cMetaDataPartiesCachees>  *  aRes = new std::vector<cMetaDataPartiesCachees>; 
   for (int aK=0 ; aK<int(mVIm->size()) ; aK++)
   {
       const std::string & aName  = (*mVIm)[aK];
       cMetaDataPartiesCachees aMDPC =
                  StdGetObjFromFile<cMetaDataPartiesCachees>
                  (
                       mWorkDir+mICNM->Assoc1To1(mCO.KeyAssocMetaData(),aName,true),
                       StdGetFileXMLSpec("SuperposImage.xml"),
                       "MetaDataPartiesCachees",
                       "MetaDataPartiesCachees"
                  );
       aRes->push_back(aMDPC);
   }

   return aRes;
}

Box2di  cAppli_Ortho::BoxImageGlob()
{
    if (mCO.BoxCalc().IsInit())
       return mCO.BoxCalc().Val();
    if (mMtDMNT) 
       return Box2di(Pt2di(0,0),mMtDMNT->NombrePixels());

   Pt2di aP0(10000000,10000000);
   Pt2di aP1 = -aP0;
   // ELISE_ASSERT(false,"cAppli_Ortho::BoxImageGlob");

   for (int aK=0 ; aK<int(mVMdpc->size()) ; aK++)
   {
       cMetaDataPartiesCachees aMDPC = (*mVMdpc)[aK];
       aP0.SetInf(aMDPC.Offset());
       aP1.SetSup(aMDPC.Offset()+aMDPC.Sz());

   }
   return Box2di(aP0,aP1);
}

const std::vector<std::string> * cAppli_Ortho::GetImNotTiled(const std::vector<std::string> * Input)
{
    std::vector<std::string> * aRes = new std::vector<std::string>;
    for (int aK=0; aK<int(Input->size()) ; aK++)
    {
        const std::string & aName = (*Input)[aK];
        if (! Tiff_Im::IsNameInternalTile(aName,mICNM))
           aRes->push_back(aName);
    }


    return aRes;
}


cAppli_Ortho::cAppli_Ortho     
(
    const cCreateOrtho & aCO,
    int argc,
    char ** argv
) :
    mWorkDir  (StdWorkdDir(aCO.WorkDir(),argv[1])),
    mCO       (aCO),
    mICNM     (cInterfChantierNameManipulateur::StdAlloc
                    (
                         argc,argv,
                         mWorkDir,
                         mCO.FileChantierNameDescripteur(),
                         mCO.DicoLoc().PtrVal()
                    )
              ),
    mVIm        (GetImNotTiled(mICNM->Get(mCO.KeySetIm()))),
    mVMdpc      (VMDPC()),
    mMtDMNT     (GetMtdMNT()),
    mTF0        (0),
    mFileOrtho  (0),
    mBoxCalc    (BoxImageGlob()),
    mImEtiqTmp  (1,1),
    mTImEtiqTmp (mImEtiqTmp),
    mImIndex    (1,1),
    mTImIndex   (mImIndex),
    mScoreNadir (1,1),
    mTScoreNadir (mScoreNadir),
    mLutLabelR  (1),
    mLutLabelV  (1),
    mLutLabelB  (1),
    mLutInd     (1),
    mImMasqMesure (1,1),
    mTImMasqMesure (mImMasqMesure),
    mW          (0),
    mEgalise    (mCO.SectionEgalisation().IsInit()),
    mSE         (mEgalise ? &(mCO.SectionEgalisation().Val()) : 0),
    mCompMesEg  (false),
    mERPrinc    (0),
    mERPhys     (0),
    mDoEgGlob   (false),
    mSeuilCorrel (-1),
    mNbPtMoyPerIm (0),
    mDynGlob      (aCO.DynGlob().Val()),
    mNbIm2Test   (0)
{
// std::cout << "DGGGGGgg "<< mDynGlob << "\n"; getchar();
    if (mEgalise)
    {
        mNameFileMesEg = mICNM->Dir()+mSE->NameFileMesures();
        mCompMesEg =     (! ELISE_fp::exist_file(mNameFileMesEg)) 
                      || (!mSE->UseFileMesure().Val())
                      || (mSE->RapelOnEgalPhys().Val()) ;
        mDoEgGlob =  mSE->GlobRappInit().DoGlob().Val();

        mSzVC = mSE->SzMaxVois().Val();
        mVoisCorrel = StdPointOfCouronnes(mSE->SzMaxVois().Val(),mSE->Use4Vois().Val());
        mSeuilCorrel = mSE->CorrelThreshold().Val();

        if (0)
        {
           for (int aKV=0 ; aKV< int(mVoisCorrel.size()) ; aKV++)
           {
               std::vector<Pt2di> aV=mVoisCorrel[aKV];
               for (int aKP=0 ; aKP< int(aV.size()) ; aKP++)
                   std::cout << aV[aKP] << " ";
                std::cout <<"\n";
           }
           getchar();
        }
    }
   
    Tiff_Im::SetDefTileFile(mCO.SzTileResult().Val());

    // const std::vector<std::string>  * aVIM = mICNM->Get(mCO.KeySetIm());

    ELISE_ASSERT(mVIm->size()!=0,"No image found for ortho");

    int aCpt = 0;
    mNbPtMoyPerIm = 0;

    for (int aK=0 ; aK<int(mVIm->size()) ; aK++)
    {
        const std::string & aName  = (*mVIm)[aK];
        cMetaDataPartiesCachees aMDPC = (*mVMdpc)[aK];

        Box2di aBox (aMDPC.Offset(),aMDPC.Offset()+aMDPC.Sz());
        if ( ! InterVide(aBox,mBoxCalc))
        {
           cOneImOrhto * anO = new cOneImOrhto(aK,*this,aName,aMDPC,aBox);

           if (anO->Im2Test()) mNbIm2Test++;

           Pt2di aSz = anO->TF().sz();
           aCpt++;
           mNbPtMoyPerIm += aSz.x * aSz.y;

           mVAllOrhtos.push_back(anO);
           if (mCO.SectionSimulImage().IsInit())
           {
              anO->DoSimul(mCO.SectionSimulImage().Val());
           }
        }
    }
    mNbPtMoyPerIm /= aCpt;
    // double aPtPerIm = 1e4;
    // double anEcart =sqrt(aNb/aPtPerIm);
    // std::cout << " NB = " << aNb << " Ecart " << anEcart << "\n";

    mLutLabelR = InitRanLutLabel();
    mLutLabelV = InitRanLutLabel();
    mLutLabelB = InitRanLutLabel();
    mLutInd = Im1D_U_INT2((int)mVAllOrhtos.size());

    mTF0 = new Tiff_Im(mVAllOrhtos[0]->TF());
    mNbCh = mTF0->nb_chan();
    for (int aK=0 ; aK<mNbCh ; aK++)
        mRadMax.push_back(0);

    bool IsNew;
    Tiff_Im aTF = Tiff_Im::CreateIfNeeded
                  (
                       IsNew,
                       mWorkDir + mCO.NameOrtho(),
                       mBoxCalc.sz(),
                       //mMtDMNT.NombrePixels(),
                       mTF0->type_el(),
                       Tiff_Im::No_Compr,
                       mTF0->phot_interp()
                  );
    mFileOrtho = new Tiff_Im(aTF);

    if (0)
    {
        Pt2di aSzW(1000,800);
        aSzW.SetInf(mBoxCalc.sz());
        std::cout << aSzW << "\n";
        mW = Video_Win::PtrWStd(aSzW);
    }


    cFileOriMnt * aFOM = GetMtdMNT();
    if (aFOM)
    {
        GenTFW(*aFOM,StdPrefix(mWorkDir + mCO.NameOrtho())+".tfw");
         //const cFileOriMnt & aFOM,const std::string & aName)

        // MTDOrtho.tfw
    }

    // inicializa cache de máscaras por imagem (vazio)
    mPerImMaskFiles.clear();
    mPerImMaskFiles.resize(mVAllOrhtos.size(), (Tiff_Im*)0);
}

Video_Win *  cAppli_Ortho::W()
{
   return mW;
}

Im1D_U_INT2 cAppli_Ortho::InitRanLutLabel()
{
   Im1D_U_INT2 aRes((int)mVAllOrhtos.size());
   ELISE_COPY(aRes.all_pts(),LabelNoIm*frandr(),aRes.out());
   return aRes;
}


void cAppli_Ortho::AddOneMasqMesure(Fonc_Num & aFMGlob,int & aCpt,const cMasqMesures & aMM,const Box2di & aBoxIn)
{
     ELISE_ASSERT( mMtDMNT!=0,"No Meta data file for Ortho created by <FileMNT>");
     Tiff_Im aTF = Tiff_Im::UnivConvStd(mWorkDir+aMM.NameFile());
     Fonc_Num aF = AdaptFonc2FileOriMnt
                          (
                              "File=" + aMM.NameFile(),
                              *mMtDMNT,
                              GetOriMnt(aMM.NameMTD()),
                              aTF.in_bool_proj(),
                              false,
                              0,
                              Pt2dr(aBoxIn._p0)
                          );
     if (aCpt==0)
        aFMGlob = aF;
     else
        aFMGlob = aFMGlob && aF;
     aCpt++;
}


void cAppli_Ortho::ResetMasqMesure(const Box2di & aBoxIn)
{
    mImMasqMesure.Resize(aBoxIn.sz());
    mTImMasqMesure =  TIm2D<U_INT1,INT>(mImMasqMesure);
    ELISE_COPY(mImMasqMesure.all_pts(),1,mImMasqMesure.out());

    Fonc_Num aFMGlob = 1;
    int aCpt=0;

    for 
    (
           std::list<cMasqMesures>::const_iterator itM=mCO.ListMasqMesures().begin();
           itM !=mCO.ListMasqMesures().end();
           itM++
    )
    {
           AddOneMasqMesure(aFMGlob,aCpt,*itM,aBoxIn);

    }
    for 
    (
           std::list<std::string>::const_iterator itS=mCO.FileExterneMasqMesures().begin();
           itS !=mCO.FileExterneMasqMesures().end();
           itS++
    )
    {
           cMasqMesures aMM = StdGetObjFromFile<cMasqMesures>
                              (
                                  mWorkDir+*itS,
                                  StdGetFileXMLSpec("SuperposImage.xml"),
                                  "MasqMesures",
                                  "MasqMesures"
                              );
           AddOneMasqMesure(aFMGlob,aCpt,aMM,aBoxIn);

    }



    if (aCpt!=0)
    {
        ELISE_COPY
        (
             mImMasqMesure.all_pts(),
             aFMGlob,
             mImMasqMesure.out()
        );
    }
}


void  cAppli_Ortho::DoOneBox
      (
          const Box2di& aBoxOut,
          const Box2di& aBoxIn,
          eModeMapBox   aMode
      )
{
    mSzCur = aBoxIn.sz();
    mCurBoxIn = aBoxIn;
    mCurBoxOut = aBoxOut;
    Resize(mSzCur);
    std::cout << "OUT " << aBoxOut._p0 <<  aBoxOut._p1 << "\n";
    std::cout << " IN " << aBoxIn._p0 <<  aBoxIn._p1 << "\n\n";

    int aKReserve =0;
    mVLI.clear();
    mNbLoadedIm2Test = 0;

    for (int aKOr = 0 ; aKOr<int(mVAllOrhtos.size()) ; aKOr++)
    {
        cOneImOrhto * anOrth = mVAllOrhtos[aKOr];
        if (anOrth->Instersecte(aBoxIn))
        {
            while (aKReserve>=int(mReserveLoadedIms.size()))
            {
               int anInd = (int)mReserveLoadedIms.size();
               mReserveLoadedIms.push_back(new cLoadedIm(*this,anInd,*mTF0,mSzMaxIn));
            }
            cLoadedIm * aLI = mReserveLoadedIms[aKReserve];

            if (aLI->Init(*anOrth,aBoxOut,aBoxIn))
            {
                mLutInd.data()[aKReserve] = anOrth->TheInd();
                aKReserve++;
                mVLI.push_back(aLI);
                if (aLI->Im2Test())
                   mNbLoadedIm2Test++;
            }
        }
    }
    ELISE_ASSERT(mVLI.size()<LabelNoIm,"Too much image in DoOneBox");

    // std::cout << "IM2TEST " <<  mNbLoadedIm2Test  << " " << mNbIm2Test << "\n";
    if (mNbLoadedIm2Test < mNbIm2Test)
    {
        std::cout <<  "Squeezed bu testing rules\n";
        return;
    }


    if (aMode==eModeCompMesSeg)
    {
        ResetMasqMesure(aBoxIn);
        DoIndexNadir();
        RemplitOneStrEgal();
    }
    if (aMode==eModeOrtho)
    {
       OrthoRedr();

       std::cout << " Radiom Max " ;
       for (int aK=0 ; aK<mNbCh ; aK++)
       {
            std::cout << " " << mRadMax[aK];
       }
       std::cout << "\n";
    }
}

/*
void cAppli_Ortho::OrthoRedr()
{
   DoIndexNadir();
   MakeOrthoOfIndex();
   VisuLabel();
   // getchar();
   SauvLabel();
   SauvAll();
}
*/

void cAppli_Ortho::OrthusCretinus()
{
   DoIndexNadir();
   MakeOrthoOfIndex();
   SauvAll();
   
}



bool cAppli_Ortho::DoEgGlob() const
{
  return mDoEgGlob;
}


void cAppli_Ortho::SauvAll()
{   
    SauvOrtho();
    SauvPerImageMasks();  // <-- salva máscaras por imagem (retângulos mBox)
    SauvLabel();
}
void cAppli_Ortho::SauvOrtho()
{
    ELISE_COPY
    (
        rectangle(mCurBoxOut._p0,mCurBoxOut._p1),
        trans(StdInput(mIms),-mCurBoxIn._p0),
        mFileOrtho->out()
    );
}

void cAppli_Ortho::SauvLabel()
{
   if ((!mCO.NameLabels().IsInit()) || (mCO.NameLabels().Val() =="NoLabel"))
      return;
   bool IsNew;
   Tiff_Im aTF = Tiff_Im::CreateIfNeeded
                  (
                       IsNew,
                       mWorkDir + mCO.NameLabels().Val(),
                       mBoxCalc.sz(),
                       GenIm::u_int1,
                       Tiff_Im::No_Compr,
                       Tiff_Im::BlackIsZero
                  );
    Symb_FNum aSL (mImIndex.in()); 

    Fonc_Num aFonc =   Min(254,(aSL>=0) * mLutInd.in()[Max(0,aSL)])
                     + 255 * (aSL<0);


    ELISE_COPY
    (
        rectangle(mCurBoxOut._p0,mCurBoxOut._p1),
        trans(aFonc,-mCurBoxIn._p0),
        aTF.out()
    );
}

void cAppli_Ortho::VisuLabel()
{
  if (!mW) return;
  ELISE_COPY(mW->all_pts(),P8COL::black,mW->odisc());
  ELISE_COPY
  (
     select(mImIndex.all_pts(),mImIndex.in()>=0),
     Virgule
     (
          mLutLabelR.in()[mImIndex.in()],
          mLutLabelV.in()[mImIndex.in()],
          mLutLabelB.in()[mImIndex.in()]
     ),
     mW->orgb()
  );
}

void cAppli_Ortho::MakeOrthoOfIndex()
{
    for (int aK=0; aK<int(mIms.size()) ; aK++)
    {
        mIms[aK]->Resize(mSzCur);
        ELISE_COPY(mIms[aK]->all_pts(),0,mIms[aK]->out());
    }
    Pt2di aP;
    for (aP.x=0 ; aP.x<mSzCur.x ; aP.x++)
    {
       for (aP.y=0 ; aP.y<mSzCur.y ; aP.y++)
       {
           int aInd = mImIndex.GetI(aP);
           if (aInd != theNoIndex)
           {
               mVLI[aInd]->TransfertIm(mIms,aP,aP+mCurBoxIn._p0);
           }
       }
    }
}
 

void cAppli_Ortho::Resize(const Pt2di & aSz)
{
   mImIndex.Resize(aSz);
   mTImIndex =  TIm2D<INT2,INT>(mImIndex);

   mScoreNadir.Resize(aSz);
   mTScoreNadir =  TIm2D<REAL4,REAL8>(mScoreNadir);

   mImEtiqTmp.Resize(aSz);
   mTImEtiqTmp = TIm2D<U_INT1,INT>(mImEtiqTmp);
   ELISE_COPY(mImEtiqTmp.all_pts(),0,mImEtiqTmp.out());
   ELISE_COPY(mImEtiqTmp.border(1),10,mImEtiqTmp.out());
}



void cAppli_Ortho::DoIndexNadir()
{

   ELISE_COPY(mScoreNadir.all_pts(),1e4,mScoreNadir.out());
   ELISE_COPY(mImIndex.all_pts(),theNoIndex,mImIndex.out());

   for (int aKI=0 ; aKI<int(mVLI.size()) ; aKI++)
   {
        mVLI[aKI]->UpdateNadirIndex(mScoreNadir,mImIndex);
   }
}


void cAppli_Ortho::DoAll()
{
    if (mEgalise)
    {
        DoEgalise();
    }
    DoOrtho();
}


void cAppli_Ortho::DoOrtho()
{
   MapBoxes(eModeOrtho);
}

void cAppli_Ortho::MapBoxes(eModeMapBox aMode)
{
   int aSz = mCO.SzDalle();
   int aBrd = mCO.SzBrd();
   if (aMode==eModeCompMesSeg)
      aBrd = 0;
   Pt2di aPBrd(aBrd,aBrd);
   cDecoupageInterv2D   aDI2D(mBoxCalc,Pt2di(aSz,aSz),Box2di(-aPBrd,aPBrd));

   mSzMaxIn = aDI2D.SzMaxIn();

   Resize(mSzMaxIn);
  
   mIms = mTF0->VecOfIm(mSzMaxIn);
   for (int aKB=mCO.KBox0().Val(); aKB<aDI2D.NbInterv() ;aKB++)
   {
       std::cout << "KBOX = " << aKB << " On " << aDI2D.NbInterv() << "\n";
       DoOneBox
       ( 
            aDI2D.KthIntervOut(aKB),
            aDI2D.KthIntervIn(aKB),
            aMode
       );
   }
   
}


cInterfChantierNameManipulateur * cAppli_Ortho::ICNM() const
{
   return mICNM;
}
const cCreateOrtho & cAppli_Ortho::CO() const
{
   return mCO;
}

const std::string &  cAppli_Ortho::WD() const
{
   return mWorkDir;
}

/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant �  la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est régi par la licence CeCILL-B soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargamento,  �  l'utilisation,  �  la modificação e/ou ao
desenvolvimento e  a reprodução do logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui pode le rendre complexe a  
manipuler e que le réserve donc a des développeurs e des profissionais
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités a charger  e  tester  l'adéquation  do
logiciel a leurs besoins dans les conditions permitindo d'assurer la
segurança de leurs systèmes et ou de leurs dados e, plus généralement, 
a l'utiliser e l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder a cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL-B, et que vous en avez accepté les
termes.
Footer-MicMac-eLiSe-25/06/2007*/


// -----------------------------------------------------------------------------
// ADDED BY GABRIEL
// Now: create per-image rectangular masks based on cOneImOrhto::mBox
// For each tile, write intersection between that image mBox and the current tile
// into the per-image TIFF mask (one TIFF per input image, size = full ortho).
// -----------------------------------------------------------------------------

Tiff_Im * cAppli_Ortho::GetPerImMaskFile(int aGlobIndex)
{
    // garante vetor do tamanho apropriado
    if ((int)mPerImMaskFiles.size() != (int)mVAllOrhtos.size())
        mPerImMaskFiles.resize(mVAllOrhtos.size(), (Tiff_Im*)0);

    if (mPerImMaskFiles[aGlobIndex])
        return mPerImMaskFiles[aGlobIndex];

    std::string dir = mWorkDir + "OrthoByFrameMasks/";
    ELISE_fp::MkDirRec(dir);
    
    // Usa nome original da imagem para gerar nome de máscara consistente
    const std::string & aOriginalName = (*mVIm)[aGlobIndex];
    std::string aMaskName = dir + "Mask_" + StdPrefix(aOriginalName) + ".tif";

    bool IsNew;
    Tiff_Im aTF = Tiff_Im::CreateIfNeeded
                  (
                      IsNew,
                      aMaskName,
                      mBoxCalc.sz(),
                      GenIm::u_int1,
                      Tiff_Im::No_Compr,
                      Tiff_Im::BlackIsZero
                  );
    
    mPerImMaskFiles[aGlobIndex] = new Tiff_Im(aTF);
    return mPerImMaskFiles[aGlobIndex];
}

void cAppli_Ortho::SauvPerImageMasks()
{
    // Se não há imagens carregadas neste tile, nada a fazer
    if (mVLI.empty()) 
        return;

    // Usaremos mImEtiqTmp (redimensionado para mSzCur em Resize) como buffer temporário.
    // Garante que está zerado
    ELISE_COPY(mImEtiqTmp.all_pts(), 0, mImEtiqTmp.out());

    int nLocal = (int)mVLI.size();

    for (int aLocal = 0; aLocal < nLocal; ++aLocal)
    {
        int aGlobIndex = mLutInd.data()[aLocal];
        if (aGlobIndex < 0 || aGlobIndex >= (int)mVAllOrhtos.size())
            continue;

        cOneImOrhto * anOrth = mVAllOrhtos[aGlobIndex];
        if (!anOrth) continue;

        // obtém bbox da imagem no espaço da ortofoto
        Box2di aImgBox = anOrth->Box(); // assume-se getter Box() existe e devolve mBox
        // interseção com a box do tile de saída (mCurBoxOut) para evitar escrita fora do tile
        Pt2di i0( std::max(aImgBox._p0.x, mCurBoxOut._p0.x), std::max(aImgBox._p0.y, mCurBoxOut._p0.y) );
        Pt2di i1( std::min(aImgBox._p1.x, mCurBoxOut._p1.x), std::min(aImgBox._p1.y, mCurBoxOut._p1.y) );

        // Se interseção vazia, nada a fazer para este tile
        if ((i0.x >= i1.x) || (i0.y >= i1.y))
            continue;

        // coords locais relativos a mCurBoxIn (buffer mImEtiqTmp)
        Pt2di loc0 = i0 - mCurBoxIn._p0;
        Pt2di loc1 = i1 - mCurBoxIn._p0;

        // Preenche a sub-rect no buffer temporário com 255
        ELISE_COPY(rectangle(loc0, loc1), 255, mImEtiqTmp.out());

        // Obtém o arquivo de máscara para a imagem e escreve a região do tile
        Tiff_Im * pMask = GetPerImMaskFile(aGlobIndex);
        if (pMask)
        {
            // Translada o buffer local para coordenadas globais e escreve no TIF da máscara
            ELISE_COPY
            (
                rectangle(mCurBoxOut._p0, mCurBoxOut._p1),
                trans(mImEtiqTmp.in(), -mCurBoxIn._p0),
                pMask->out()
            );
        }

        // limpa a sub-rect no buffer temporário para próxima imagem
        ELISE_COPY(rectangle(loc0, loc1), 0, mImEtiqTmp.out());
    }
}

cAppli_Ortho::~cAppli_Ortho()
{
    // libera eventuais Tiff_Im* alocados para máscaras
    for (size_t i = 0; i < mPerImMaskFiles.size(); ++i)
    {
        if (mPerImMaskFiles[i] != 0)
        {
            delete mPerImMaskFiles[i];
            mPerImMaskFiles[i] = 0;
        }
    }
}
// END ADDED BY GABRIEL
