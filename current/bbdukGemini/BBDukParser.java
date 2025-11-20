package bbdukGemini;

import java.io.File;
import java.util.ArrayList;
import shared.Parse;
import shared.Tools;

/**
 * Handles parsing for BBDuk-specific arguments.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 19, 2025
 */
public class BBDukParser{
	
	public BBDukParser(){}
	
	public boolean parse(String arg, String a, String b, int i){
		
		if(a.equals("outb") || a.equals("outm") || a.equals("outb1") || a.equals("outm1") ||
		   a.equals("outbad") || a.equals("outbad1") || a.equals("outmatch") || a.equals("outmatch1")){
			outb1=b;
		}else if(a.equals("outb2") || a.equals("outm2") || a.equals("outbad2") || a.equals("outmatch2")){
			outb2=b;
		}else if(a.equals("qfoutb") || a.equals("qfoutm") || a.equals("qfoutb1") || a.equals("qfoutm1")){
			qfoutb1=b;
		}else if(a.equals("qfoutb2") || a.equals("qfoutm2")){
			qfoutb2=b;
		}
		
		else if(a.equals("ottm") || a.equals("outputtrimmedtomatch")){
			addTrimmedToBad=Parse.parseBoolean(b);
		}else if(a.equals("rename")){
			rename=Parse.parseBoolean(b);
		}else if(a.equals("refnames") || a.equals("userefnames")){
			useRefNames=Parse.parseBoolean(b);
		}else if(a.equals("ordered")){
			ordered=Parse.parseBoolean(b);
		}
		
		else if(a.equals("ref") || a.equals("adapters")){
			ref=(b==null) ? null : (new File(b).exists() ? new String[]{b} : b.split(","));
		}else if(a.equals("altref")){
			altref=(b==null) ? null : (new File(b).exists() ? new String[]{b} : b.split(","));
		}else if(a.equals("literal")){
			literal=BBDukParser.processLiteralArg(b);
		}else if(a.equals("k")){
			assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=27\n";
			k=Integer.parseInt(b);
		}else if(a.equals("mink") || a.equals("kmin")){
			mink=Integer.parseInt(b);
		}else if(a.equals("useshortkmers") || a.equals("shortkmers") || a.equals("usk")){
			useShortKmers=Parse.parseBoolean(b);
		}else if(a.equals("trimextra") || a.equals("trimpad") || a.equals("tp")){
			trimPad=Integer.parseInt(b);
		}else if(a.equals("hdist") || a.equals("hammingdistance")){
			hammingDistance=Integer.parseInt(b);
		}else if(a.equals("qhdist") || a.equals("queryhammingdistance")){
			qHammingDistance=Integer.parseInt(b);
		}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
			editDistance=Integer.parseInt(b);
		}else if(a.equals("hdist2") || a.equals("hammingdistance2")){
			hammingDistance2=Integer.parseInt(b);
		}else if(a.equals("qhdist2") || a.equals("queryhammingdistance2")){
			qHammingDistance2=Integer.parseInt(b);
		}else if(a.equals("edits2") || a.equals("edist2") || a.equals("editdistance2")){
			editDistance2=Integer.parseInt(b);
		}else if(a.equals("maxskip") || a.equals("maxrskip") || a.equals("mxs")){
			maxSkip=Integer.parseInt(b);
		}else if(a.equals("minskip") || a.equals("minrskip") || a.equals("mns")){
			minSkip=Integer.parseInt(b);
		}else if(a.equals("skip") || a.equals("refskip") || a.equals("rskip")){
			minSkip=maxSkip=Integer.parseInt(b);
		}else if(a.equals("mm") || a.equals("maskmiddle")){
			if(b==null || Tools.startsWithLetter(b)){
				maskMiddle=Parse.parseBoolean(b);
			}else{
				midMaskLen=Integer.parseInt(b);
				maskMiddle=midMaskLen>0;
			}
		}else if(a.equals("rcomp")){
			rcomp=Parse.parseBoolean(b);
		}else if(a.equals("forbidns") || a.equals("forbidn") || a.equals("fn")){
			forbidNs=Parse.parseBoolean(b);
		}else if(a.equals("findbestmatch") || a.equals("fbm")){
			findBestMatch=Parse.parseBoolean(b);
		}
		
		else if(a.equals("ktrim")){
			if(b==null){b="";}
			if(b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr") || b.equalsIgnoreCase("tips")){
				ktrimLeft=ktrimRight=true; ktrimN=ksplit=false;
			}else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){
				ktrimLeft=true; ktrimRight=false; ktrimN=ksplit=false;
			}else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){
				ktrimLeft=false; ktrimRight=true; ktrimN=ksplit=false;
			}else if(b.equalsIgnoreCase("n")){
				ktrimLeft=ktrimRight=ksplit=false; ktrimN=true;
			}else{
				ktrimRight=ktrimLeft=ktrimN=false;
			}
		}else if(a.equals("ksplit")){
			ksplit=Parse.parseBoolean(b);
			if(ksplit){ktrimLeft=ktrimRight=ktrimN=false;}
		}else if(a.equals("ktrimexclusive")){
			ktrimExclusive=Parse.parseBoolean(b);
		}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
			trimByOverlap=Parse.parseBoolean(b);
		}else if(a.equals("strictoverlap")){
			strictOverlap=Parse.parseBoolean(b);
		}else if(a.equals("minoverlap")){
			minOverlap=Integer.parseInt(b);
		}else if(a.equals("mininsert")){
			minInsert=Integer.parseInt(b);
		}
		
		else if(a.equals("mingc")){
			minGC=Float.parseFloat(b);
		}else if(a.equals("maxgc")){
			maxGC=Float.parseFloat(b);
		}else if(a.equals("filtergc")){
			filterGC=Parse.parseBoolean(b);
		}
		
		else if(a.equals("minbasefrequency")){
			minBaseFrequency=Float.parseFloat(b);
		}else if(a.equals("mincoveredfraction") || a.equals("mincovfraction") || a.equals("mcf")){
			minCoveredFraction=Float.parseFloat(b);
		}else if(a.equals("minkmerfraction") || a.equals("minfraction") || a.equals("mkf")){
			minKmerFraction=Float.parseFloat(b);
		}else if(a.equals("maxbadkmers") || a.equals("mbk")){
			maxBadKmers=Integer.parseInt(b);
		}else if(a.equals("pratio") || a.equals("polymerratio")){
			countPolymers=true;
		}else if(a.equals("entropy") || a.equals("minentropy")){
			entropyCutoff=Float.parseFloat(b);
			calcEntropy=(entropyCutoff>=0);
		}else if(a.equals("entropymask")){
			entropyMask=Parse.parseBoolean(b);
			if(entropyMask){calcEntropy=true;}
		}else if(a.equals("entropytrim")){
			entropyTrim=Parse.parseBoolean(b) ? 3 : 0;
			if(entropyTrim>0){calcEntropy=true;}
		}else if(a.equals("filtervars")){
			filterVars=Parse.parseBoolean(b);
		}else if(a.equals("minx")){ xMinLoc=Integer.parseInt(b); locationFilter=true;
		}else if(a.equals("maxx")){ xMaxLoc=Integer.parseInt(b); locationFilter=true;
		}else if(a.equals("miny")){ yMinLoc=Integer.parseInt(b); locationFilter=true;
		}else if(a.equals("maxy")){ yMaxLoc=Integer.parseInt(b); locationFilter=true;
		}else if(a.equals("maxbasesoutu")){ maxBasesOutu=Parse.parseKMG(b);
		}else if(a.equals("maxbasesoutm")){ maxBasesOutm=Parse.parseKMG(b);
		}else if(a.equals("pairedtosingle")){ pairedToSingle=Parse.parseBoolean(b);
		}else if(a.equals("trimfailures") || a.equals("trimfailuresto1bp")){
			trimFailuresTo1bp=Parse.parseBoolean(b);
		}else if(a.equals("restrictleft")){
			restrictLeft=Integer.parseInt(b);
		}else if(a.equals("restrictright")){
			restrictRight=Integer.parseInt(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public void postParse(boolean amino){
		if(hammingDistance2==-1){hammingDistance2=hammingDistance;}
		if(qHammingDistance2==-1){qHammingDistance2=qHammingDistance;}
		if(editDistance2==-1){editDistance2=editDistance;}
		hammingDistance=Tools.max(editDistance, hammingDistance);
		hammingDistance2=Tools.max(editDistance2, hammingDistance2);
		
		final int maxSupportedK=(amino ? 12 : 31);
		kbig=(k>maxSupportedK ? k : -1);
		k=Tools.min(k, maxSupportedK);
		k2=k-1; 
		
		if(mink<1){mink=6;}
		mink=Tools.min(mink, k);
		
		kfilter=(ref!=null || literal!=null) && !(ktrimRight || ktrimLeft || ktrimN || ksplit);
		if(mink>0 && mink<k){useShortKmers=true;}
		
		minSkip=Tools.max(1, Tools.min(minSkip, maxSkip));
		maxSkip=Tools.max(minSkip, maxSkip);
		
		if(maskMiddle){
			midMaskLen=(midMaskLen>0 ? midMaskLen : 2-(k&1));
		}else{
			midMaskLen=0;
		}
		
		if(kbig>k){
			minSkip=maxSkip=0;
			maskMiddle=false;
			midMaskLen=0;
			useShortKmers=false;
		}
		
		if((ktrimLeft || ktrimRight || ktrimN || ksplit) && kbig>k){
			kbig=k;
		}
		
		if(strictOverlap){
			maxRatio=0.05f;
			ratioMargin=9f;
			ratioOffset=0.5f;
			efilterRatio=3.5f;
			efilterOffset=0.05f;
			pfilterRatio=0.001f;
			meeFilter=15f;
		}else{
			maxRatio=0.10f;
			ratioMargin=5f;
			ratioOffset=0.4f;
			efilterRatio=6f;
			efilterOffset=0.05f;
			pfilterRatio=0.00005f;
			meeFilter=999999999;
		}

		if(minOverlap<0){minOverlap=minOverlap0;}
		if(minInsert<0){minInsert=minInsert0;}
		
		if(minGC>0 || maxGC<1){filterGC=true;}

		BBDukIndexConstants.initialize(k, mink, amino, midMaskLen);
	}
	
	// ... [static methods omitted for brevity] ...
	public static final String[] processLiteralArg(String arg){
		if(arg==null){return null;}
		String[] split=arg.split(",");
		ArrayList<String> list=new ArrayList<String>(split.length);
		for(String b : split){
			String c=processLiteralTerm(b);
			if(c!=null){list.add(c);}
		}
		return list.isEmpty() ? null : list.toArray(new String[0]);
	}
	
	public static final String processLiteralTerm(String b){
		assert(b.length()>0) : "Invalid literal sequence: '"+b+"'";
		if(dna.AminoAcid.isACGTN(b)){return b;}
		b=b.toUpperCase().replaceAll("-", "");
		if(b.startsWith("POLY")){
			b=b.replace("POLY", "");
			assert(dna.AminoAcid.isACGTN(b)) : "Invalid literal sequence: '"+b+"'";
			StringBuilder sb=new StringBuilder(40);
			final int minlen=Tools.max(31+b.length(), b.length()*3);
			while(sb.length()<minlen){sb.append(b);}
			return sb.toString();
		}
		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public String[] ref=null;
	public String[] altref=null;
	public String[] literal=null;
	
	public String outb1=null;
	public String outb2=null;
	public String qfoutb1=null;
	public String qfoutb2=null;
	
	public boolean addTrimmedToBad=false;
	public boolean rename=false;
	public boolean useRefNames=false;
	public boolean histogramsBeforeProcessing=true;
	public boolean ordered=false;
	public long maxBasesOutu=-1;
	public long maxBasesOutm=-1;
	
	public int k=31;
	public int kbig=-1;
	public int k2=30;
	public int mink=-1;
	public boolean useShortKmers=false;
	public boolean maskMiddle=true;
	public int midMaskLen=0;
	public boolean rcomp=true;
	public boolean forbidNs=false;
	public float minKmerFraction=0;
	public int maxBadKmers=0;
	
	public int hammingDistance=0;
	public int qHammingDistance=0;
	public int editDistance=0;
	public int hammingDistance2=-1;
	public int qHammingDistance2=-1;
	public int editDistance2=-1;
	
	public int maxSkip=1;
	public int minSkip=1;
	public int trimPad=0;
	public byte trimSymbol='N';
	public boolean kmaskFullyCovered=false;
	public boolean kmaskLowercase=false;
	
	public boolean kfilter=false;
	public boolean ktrimLeft=false;
	public boolean ktrimRight=false;
	public boolean ktrimN=false;
	public boolean ksplit=false;
	public boolean ktrimExclusive=false;
	public boolean findBestMatch=false;
	
	public float minCoveredFraction;
	public float minBaseFrequency;
	
	public int forceTrimModulo=-1;
	public int forceTrimLeft=-1;
	public int forceTrimRight=-1;
	public int forceTrimRight2=-1;
	public boolean tossJunk=false;
	public boolean trimFailuresTo1bp=false;
	public boolean removePairsIfEitherBad=true;
	public boolean pairedToSingle=false;
	public boolean filterVars=false;
	
	public boolean trimByOverlap=false;
	public boolean strictOverlap=false;
	public boolean useQualityForOverlap=false;
	public int minOverlap0=7;
	public int minOverlap=-1;
	public int minInsert0=16;
	public int minInsert=-1;
	public float maxRatio;
	public float ratioMargin;
	public float ratioOffset;
	public float efilterRatio;
	public float efilterOffset;
	public float pfilterRatio;
	public float meeFilter;
	
	public boolean swift=false;
	public int trimPolyA=0;
	public int trimPolyGLeft=0;
	public int trimPolyGRight=0;
	public int trimPolyCLeft=0;
	public int trimPolyCRight=0;
	public int filterPolyG=0;
	public int filterPolyC=0;
	
	public boolean calcEntropy=false;
	public boolean entropyMask=false;
	public int entropyTrim=0;
	public boolean entropyMark=false;
	public float entropyCutoff=-1;
	public boolean entropyHighpass=true;
	public boolean entropyMaskLowercase=false;
	
	public boolean countPolymers=false;
	public boolean quantizeQuality=false;
	
	public boolean locationFilter=false;
	public int xMinLoc=-1, xMaxLoc=-1, yMinLoc=-1, yMaxLoc=-1;
	public int restrictLeft=0, restrictRight=0;
	
	public float minGC=0;
	public float maxGC=1;
	public boolean filterGC=false;
	public int maxBadSubs=0;

}