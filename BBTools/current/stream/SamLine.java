package stream;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import dna.ScafLoc;
import shared.KillSwitch;
import shared.LineParser1;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import var2.ScafMap;
import var2.Scaffold;


public class SamLine implements Serializable {
	
//	426_647_582	161	chr1	10159	0	26M9H	chr3	170711991	0	TCCCTAACCCTAACCCTAACCTAACC	IIFIIIIIIIIIIIIIIIIIICH2<>	RG:Z:20110708003021394	NH:i:3	CM:i:2	SM:i:1	CQ:Z:A9?(BB?:<A?>=>B67=:7A);.%8'%))/%*%'	CS:Z:G12002301002301002301023010200000003	XS:A:+
	
//	1 QNAME String [!-?A-~]f1,255g Query template NAME
//	2 FLAG Int [0,216-1] bitwise FLAG
//	3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
//	4 POS Int [0,229-1] 1-based leftmost mapping POSition
//	5 MAPQ Int [0,28-1] MAPping Quality
//	6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
//	7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next fragment
//	8 PNEXT Int [0,229-1] Position of the mate/next fragment
//	9 TLEN Int [-229+1,229-1] observed Template LENgth
//	10 SEQ String \*|[A-Za-z=.]+ fragment SEQuence
//	11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
	
	
//	FCB062MABXX:1:1101:1177:2115#GGCTACAA	147	chr11	47765857	29	90M	=	47765579	-368	CCTCTGTGGCCCGGGTTGGAGTGCAGTGTCATGATCATGGCTCGCTGTAGCTACACCCTTCTGAGCTCAAGCAATCCTCCCACCTCTCCC	############################################################A@@><D<AAAB<=A2BD/BC<7:<4<%679	XT:A:M	NM:i:5	SM:i:29	AM:i:29	XM:i:5	XO:i:0	XG:i:0	MD:Z:7T4A15G26A30A3
//	FCB062MABXX:1:1101:1193:2122#GGCTACAA	77	*	    0	         0	*	*	0	           0	TATATATGTGCTATGTACAGCATTGGAATTCACACCCTACACTTTCAAAAGNGAGCCCTAAATAAATGTTAGATCGGAAGAGCACACGTC	FCFCFDDDADDEDEBDAEDFEDEFFGGFGGHEEFHHHHHHEDDDEDFFEFB#CBBA@B8BGGFGEEEC>DGGGDFBGGGGHHHHH9<@##
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -4180486051387471116L;

	public SamLine() {}
	
	public SamLine(String s){//TODO: Change to LineParser version
		this(s.split("\t"));
	}

	public SamLine(SamLine sl){
		setFrom(sl);
	}
	
//	/** Prevents references to original string, in case of e.g. very long MD tags. */
//	public SamLine toSamLine(String s){
//		String[] split=s.split("\t");
//		split[0]=new String(split[0]);
//		split[5]=new String(split[5]);
//		split[9]=new String(split[9]);
//		split[10]=new String(split[10]);
//		for(int i=11; i<split.length; i++){
//			split[i]=new String(split[i]);
//		}
//		return new SamLine(split);
//	}
	
	private void setFrom(SamLine sl){
		qname=sl.qname;
		flag=sl.flag;
		rname=sl.rname;
		rnameS=sl.rnameS;
		pos=sl.pos;
		mapq=sl.mapq;
		cigar=sl.cigar;
		rnext=sl.rnext;
		pnext=sl.pnext;
		tlen=sl.tlen;
		seq=sl.seq;
		qual=sl.qual;
		optional=sl.optional;
	}
	
	
	public SamLine(Read r1, int fragNum){
		
		if(verbose){
			System.err.println("new SamLine for read with match "+(r1.match==null ? "null" : new String(r1.match)));
		}
		
		Read r2=r1.mate;
		final boolean perfect=r1.perfect();
		
		if(Data.scaffoldLocs==null && r1.samline!=null){
			assert(SET_FROM_OK) : "Sam format cannot be used as input to this program when no genome build is loaded.\n" +
					"Please index the reference first and rerun with e.g. 'build=1', or use a different input format.";
			setFrom(r1.samline);
			return;
		}
		
//		qname=r.id.replace(' ', '_').replace('\t', '_');
//		qname=r.id.split("\\s+")[0];
		qname=r1.id.replace('\t', '_');
//		if(!KEEP_NAMES && qname.length()>2 && r2!=null){
//			if(qname.endsWith("/1") || qname.endsWith("/2") || qname.endsWith(" 1") || qname.endsWith(" 2")){}
//		}
		
		if(!KEEP_NAMES && qname.length()>2 && r2!=null){
			char c=qname.charAt(qname.length()-2);
			int num=(qname.charAt(qname.length()-1))-'1';
			if((num==0 || num==1) && (c==' ' || c=='/')){qname=qname.substring(0, qname.length()-2);}
//			if(r.pairnum()==num && (c==' ' || c=='/')){qname=qname.substring(0, qname.length()-2);}
		}
//		flag=Integer.parseInt(s[1]);
		
		int idx1=-1, idx2=-1;
		int chrom1=-1, chrom2=-1;
		int start1=-1, start2=-1, a1=0, a2=0;
		int stop1=-1, stop2=-1, b1=0, b2=0;
		int scaflen=0, scafloc=0, scaflen2=0;
		byte[] name1=bytestar, name2=bytestar;
		if(r1.mapped()){
			assert(r1.chrom>=0) : r1.chrom+", "+r1.start+", "+r1.stop;
			chrom1=r1.chrom;
			start1=r1.start;
			stop1=r1.stop;
			if(Data.isSingleScaffold(chrom1, start1, stop1)){
				assert(Data.scaffoldLocs!=null) : "\n\n"+r1+"\n\n"+r1.obj+"\n\n";
				idx1=Data.scaffoldIndex(chrom1, (start1+stop1)/2);
				name1=Data.scaffoldNames[chrom1][idx1];
				scaflen=Data.scaffoldLengths[chrom1][idx1];
				scafloc=Data.scaffoldLocs[chrom1][idx1];
				a1=Data.scaffoldRelativeLoc(chrom1, start1, idx1);
				b1=a1-start1+stop1;
			}else{
				if(verbose){System.err.println("------------- Found multi-scaffold alignment! -------------");}
				r1.setMapped(false);
				r1.setPaired(false);
				r1.match=null;
				if(r2!=null){r2.setPaired(false);}
			}
		}
		if(r2!=null && r2.mapped()){
			chrom2=r2.chrom;
			start2=r2.start;
			stop2=r2.stop;
			if(Data.isSingleScaffold(chrom2, start2, stop2)){
				idx2=Data.scaffoldIndex(chrom2, (start2+stop2)/2);
				name2=Data.scaffoldNames[chrom2][idx2];
				scaflen2=Data.scaffoldLengths[chrom2][idx2];
				a2=Data.scaffoldRelativeLoc(chrom2, start2, idx2);
				b2=a2-start2+stop2;
			}else{
				if(verbose){System.err.println("------------- Found multi-scaffold alignment for r2! -------------");}
				r2.setMapped(false);
				r2.setPaired(false);
				r2.match=null;
				if(r1!=null){r1.setPaired(false);}
			}
		}
		
		final boolean sameScaf=(r2!=null && idx1>-1 && idx1==idx2 && r1.chrom==r2.chrom);
		flag=makeFlag(r1, r2, fragNum, sameScaf);
		
		rname=r1.mapped() ? name1 : ((r2!=null && r2.mapped()) ? name2 : null);
		
		{
			int pos0, pos0_mate; //start pos
			int pos1, pos1_mate; //stop pos
			
			if(r1.mapped()){
//				int leadingClip=countLeadingClip(cigar);
				int clip=countLeadingClip(r1.match);
				int clippedIndels=countLeadingIndels(a1, r1.match);
				int tclip=countTrailingClip(r1.match);
				int tclippedIndels=countTrailingIndels(b1, scaflen, r1.match);
				
				if(verbose){
					System.err.println("leadingClip="+clip);
					System.err.println("clippedDels="+clippedIndels);
				}
				pos0=(a1+1)+clip+clippedIndels;
				pos1=(b1+1)-tclip-tclippedIndels;
				if(pos1>scaflen){pos1=scaflen;}
				
				if(pos0<1){
					//This is necessary to prevent mapped reads from having POS less than 1.
					pos0=1;
				}
				assert(pos1>=pos0) : pos0+", "+pos1+"\n"+r1+"\n"+r2+"\n";
				
			}else{
				pos0=0;
				pos1=0;
			}
			
			if(r2!=null && r2.mapped()){
				int clip=countLeadingClip(r2.match);
				int clippedIndels=countLeadingIndels(a2, r2.match);
				int tclip=countTrailingClip(r2.match);
				int tclippedIndels=countTrailingIndels(b2, scaflen, r2.match);
				if(verbose){
					System.err.println("leadingClip="+clip);
					System.err.println("clippedDels="+clippedIndels);
				}
				pos0_mate=(a2+1)+clip+clippedIndels;
				pos1_mate=(b2+1)-tclip-tclippedIndels;
				if(pos1_mate>scaflen){pos1=scaflen;}
				
				if(pos0_mate<1){
					//This is necessary to prevent mapped reads from having POS less than 1.
					pos0_mate=1;
				}
				assert(!sameScaf || pos1_mate>=pos0_mate) : pos0_mate+", "+pos1_mate+", "+scaflen+"\n"+r1+"\n"+r2+"\n";
				
			}else{
				pos0_mate=0;
				pos1_mate=0;
			}
			
			if(r2==null){
				pos=pos0;
				pnext=pos0_mate;
				tlen=0;
				assert(((pos>0 && r1.mapped()) || (pos==0 && !r1.mapped())) && pnext==0);
			}else{
				if(r1.mapped() && r2.mapped()){
					pos=pos0;
					pnext=pos0_mate;
					if(sameScaf){
//						tlen=1+(Data.max(r.stop, r2.stop)-Data.min(r.start, r2.start));
						tlen=1+(Data.max(pos1, pos1_mate)-Data.min(pos0, pos0_mate));
					}else{
						tlen=0;
					}
					assert(pos>0) : pos+"\n"+r1+"\n"+r2;
					assert(pnext>0) : pnext+"\n"+r1+"\n"+r2;
				}else if(r1.mapped() && !r2.mapped()){
					pos=pos0;
					pnext=pos0;
					tlen=0;
					assert(pos>0 && pnext>0);
				}else if(!r1.mapped() && r2.mapped()){
					pos=pos0_mate;
					pnext=pos0_mate;
					tlen=0;
					assert(pos>0 && pnext>0);
				}else if(!r1.mapped() && !r2.mapped()){
					pos=pos0;
					pnext=pos0_mate;
					tlen=0;
					assert(pos==0 && pnext==0);
				}else{assert(false);}
			}
			
			assert(pos>=0) : "Negative coordinate "+pos+" for read:\n\n"+r1+"\n\n"+r2+"\n\n"+this+"\n\na1="+a1+", a2="+a2+
				", pos0="+pos0+", pos0_mate="+pos0_mate+", clip="+countLeadingClip(cigar, true, false)+", clipM="+countLeadingClip(r1.match);
			assert(pnext>=0) : "Negative coordinate "+pnext+" for mate:\n\n"+r1+"\n\n"+r2+"\n\n"+this+"\n\na1="+a1+", a2="+a2+
				", pos0="+pos0+", pos0_mate="+pos0_mate+", clip="+countLeadingClip(cigar, true, false);
		}
		
		mapq=toMapq(r1, null);

		if(verbose){
			System.err.println("Making cigar for "+(r1.match==null ? "null" : new String(r1.match)));
		}

		final boolean inbounds=!r1.mapped() ? false : (a1>=0 && b1<scaflen);
		final boolean inbounds2=(r2==null ? true : !r2.mapped() ? false : (a2>=0 && b2<scaflen2));
		if(r1.bases!=null && r1.mapped() && r1.match!=null){
			if(VERSION>1.3f){
				if(inbounds && perfect && !r1.containsNonM()){//r.containsNonM() should be unnecessary...  it's there in case of clipping...
					cigar=(r1.length()+"=");
//					System.err.println("SETTING cigar14="+cigar);
//
//					byte[] match=r.match;
//					if(r.shortmatch()){match=Read.toLongMatchString(match);}
//					cigar=toCigar13(match, a1, b1, scaflen, r.bases);
//					System.err.println("RESETTING cigar14="+cigar+" from toCigar14("+new String(Read.toShortMatchString(match))+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}else{
					byte[] match=r1.match;
					if(r1.shortmatch()){match=Read.toLongMatchString(match);}
					cigar=toCigar14(match, a1, b1, scaflen, r1.bases);
//					System.err.println("CALLING toCigar14("+Read.toShortMatchString(match)+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}
			}else{
				if(inbounds && (perfect || !r1.containsNonNMS())){
					cigar=(r1.length()+"M");
//					System.err.println("SETTING cigar13="+cigar);
//
//					byte[] match=r.match;
//					if(r.shortmatch()){match=Read.toLongMatchString(match);}
//					cigar=toCigar13(match, a1, b1, scaflen, r.bases);
//					System.err.println("RESETTING cigar13="+cigar+" from toCigar13("+new String(Read.toShortMatchString(match))+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}else{
					byte[] match=r1.match;
					if(r1.shortmatch()){match=Read.toLongMatchString(match);}
					cigar=toCigar13(match, a1, b1, scaflen, r1.bases);
//					System.err.println("CALLING toCigar13("+Read.toShortMatchString(match)+", "+a1+", "+b1+", "+scaflen+", "+r.bases+")");
				}
			}
		}
		
		if(verbose){
			System.err.println("cigar="+cigar);
		}
		
//		assert(false);
		
//		assert(primary() || cigar.equals(stringstar)) : cigar;
//		if(pos<0){pos=0;cigar=null;rname=bytestar;mapq=0;flag|=0x4;}
		
//		assert(false) : "\npos="+pos+"\ncigar='"+cigar+"'\nVERSION="+VERSION+"\na1="+a1+", b1="+b1+"\n\n"+r.toString();
		
//		rnext=(r2==null ? stringstar : (r.mapped() && !r2.mapped()) ? "chr"+Gene.chromCodes[r.chrom] : "chr"+Gene.chromCodes[r2.chrom]);
		rnext=((r2==null || (!r1.mapped() && !r2.mapped())) ? bytestar : (r1.mapped() && r2.mapped()) ? (sameScaf ? byteequals : name2) : byteequals);
		
		assert(rnext!=byteequals || name1==name2 || name1==bytestar || name2==bytestar) :
			new String(rname)+", "+new String(rnext)+", "+new String(name1)+", "+new String(name2)+"\n"+r1+"\n"+r2;
		
//		assert(r1.pairnum()==0) : r1.mapped()+", "+r2.mapped()+"fragNum="+fragNum+
//			"\nname1="+new String(name1)+"\nname2="+new String(name2)+"\nrname="+new String(rname)+"\nrnext="+new String(rnext)+
//			"\nname1="+name1+"\nname2="+name2+"\nrname="+rname+"\nrnext="+rnext+"\nidx1="+idx1+"\nidx2="+idx2;
		
		if(Data.scaffoldPrefixes){
			 if(rname!=null && rname!=bytestar){
				 int k=Tools.indexOf(rname, (byte)'$');
				 rname=KillSwitch.copyOfRange(rname, k+1, rname.length);
			 }
			 if(rnext!=null && rnext!=bytestar){
				 int k=Tools.indexOf(rnext, (byte)'$');
				 rnext=KillSwitch.copyOfRange(rnext, k+1, rnext.length);
			 }
		}
		
//		if(r2==null || r.stop<=r2.start){
//			//plus sign
//		}else if(r2.stop<=r.start){
//			//minus sign
//			tlen=-tlen;
//		}else{
//			//They overlap... a lot.  Physically shorter than read length.
//			if(r.start<=r2.start){
//
//			}else{
//				tlen=-tlen;
//			}
//		}
		//This version is less technically correct (does not account for very short insert reads) but probably more in line with what is expected
		if(r2==null || r1.start<r2.start || (r1.start==r2.start && r1.pairnum()==0)){
			//plus sign
		}else{
			//minus sign
			tlen=-tlen;
		}
		
//		if(r.secondary()){
////			seq=qual=stringstar;
//			seq=qual=bytestar;
//		}else{
//			if(r.strand()==Gene.PLUS){
////				seq=new String(r.bases);
//				seq=r.bases.clone();
//				if(r.quality==null){
////					qual=stringstar;
//					qual=bytestar;
//				}else{
////					StringBuilder q=new StringBuilder(r.quality.length);
////					for(byte b : r.quality){
////						q.append((char)(b+33));
////					}
////					qual=q.toString();
//					qual=new byte[r.quality.length];
//					for(int i=0, j=qual.length-1; i<qual.length; i++, j--){
//						qual[i]=(byte)(r.quality[j]+33);
//					}
//				}
//			}else{
////				seq=new String(AminoAcid.reverseComplementBases(r.bases));
//				seq=AminoAcid.reverseComplementBases(r.bases);
//				if(r.quality==null){
////					qual=stringstar;
//					qual=bytestar;
//				}else{
////					StringBuilder q=new StringBuilder(r.quality.length);
////					for(int i=r.quality.length-1; i>=0; i--){
////						q.append((char)(r.quality[i]+33));
////					}
////					qual=q.toString();
//					qual=new byte[r.quality.length];
//					for(int i=0, j=qual.length-1; i<qual.length; i++, j--){
//						qual[i]=(byte)(r.quality[j]+33);
//					}
//				}
//			}
//		}
		
		if(r1.secondary() && SECONDARY_ALIGNMENT_ASTERISKS){
//			seq=qual=bytestar;
			seq=qual=null;
		}else{
			seq=r1.bases;
			if(r1.quality==null){
//				qual=bytestar;
				qual=null;
			}else{
				qual=r1.quality;
			}
		}
		
		trimNames();
		optional=makeOptionalTags(r1, r2, perfect, scafloc, scaflen, inbounds, inbounds2);
//		assert(r.pairnum()==1) : "\n"+r.toText(false)+"\n"+this+"\n"+r2;
	}
	
	public SamLine(String[] s){
		assert(!s[0].startsWith("@")) : "Tried to make a SamLine from a header: "+s[0];
		assert(s.length>=11) : "\nNot all required fields are present: "+s.length+"\nline='"+Arrays.toString(s)+"'\n";
		if(s.length<11){
			System.err.println("Invalid SamLine: "+Arrays.toString(s));
			return;
		}
		qname=s[0];
		flag=Integer.parseInt(s[1]);
		if(RNAME_AS_BYTES){
			rname=s[2].getBytes();
		}else{
			rnameS=s[2];
		}
		pos=Integer.parseInt(s[3]);
//		try {
//			Integer.parseInt(s[4]);
//		} catch (NumberFormatException e) {
//			System.err.println(Arrays.toString(s));
//		}
		mapq=Tools.isDigit(s[4].charAt(0)) ? Integer.parseInt(s[4]) : 99; //Added for non-compliant mappers that put * here
		cigar=s[5];
		rnext=s[6].getBytes();
		pnext=(s[7].charAt(0)=='*' ? 0 : Integer.parseInt(s[7]));
		tlen=Tools.isDigit(s[8].charAt(0)) ? Integer.parseInt(s[8]) : 0; //Added for non-compliant mappers that put * here
//		seq=s[9];
//		qual=s[10];
		seq=(s[9].equals(stringstar) ? null : s[9].getBytes());
		qual=(s[10].equals(stringstar) ? null : s[10].getBytes());
		
		if(mapped() && strand()==Shared.MINUS){
			if(seq!=bytestar){AminoAcid.reverseComplementBasesInPlace(seq);}
			if(qual!=bytestar){Tools.reverseInPlace(qual);}
		}
		
		if(qual!=null && qual!=bytestar){
			for(int i=0; i<qual.length; i++){qual[i]-=33;}
		}
		
		if(!PARSE_OPTIONAL || s.length<=11){return;}
		
		if(PARSE_OPTIONAL_MD_ONLY){
			optional=new ArrayList<String>(1);
			for(int i=11; i<s.length; i++){
				if(s[i].startsWith("MD:")){
					mdTag=s[i].substring(3).getBytes();
					optional.add(s[i]);
					break;
				}
			}
		}else{
			optional=new ArrayList<String>(s.length-11);
			for(int i=11; i<s.length; i++){
				if(s[i].startsWith("MD:")){
					mdTag=s[i].substring(3).getBytes();
				}
				optional.add(s[i]);
			}
		}
		
		trimNames();
	}
	
	public SamLine(LineParser1 lp){
		assert(!lp.startsWith('@')) : "Tried to make a SamLine from a header: "+lp.toString();
		
		if(PARSE_0){qname=lp.parseString(0);}
		flag=lp.parseInt(1);
		if(PARSE_2){
			if(RNAME_AS_BYTES){
				rname=(lp.termEquals('*', 2)) ? null : lp.parseByteArray(2);
			}else {
				rnameS=(lp.termEquals('*', 2)) ? null : lp.parseString(2);
			}
		}
		pos=lp.parseInt(3);
		mapq=lp.parseInt(4);
		if(PARSE_5){cigar=lp.parseString(5);}
		if(PARSE_6){rnext=lp.parseByteArray(6);}
		if(PARSE_7){pnext=lp.parseInt(7);}
		if(PARSE_8){tlen=lp.parseInt(8);}
		seq=lp.parseByteArray(9);
		if(PARSE_10){qual=lp.parseByteArray(10);}

		assert((seq==bytestar)==(Tools.equals(seq, bytestar)));
		assert((qual==bytestar)==(Tools.equals(qual, bytestar)));
		
		if(mapped() && strand()==Shared.MINUS && FLIP_ON_LOAD){
			if(seq!=bytestar){AminoAcid.reverseComplementBasesInPlace(seq);}
			if(qual!=bytestar){Tools.reverseInPlace(qual);}
		}
		
		if(qual!=null && qual!=bytestar){
			for(int i=0; i<qual.length; i++){qual[i]-=33;}
		}
		
		if(PARSE_OPTIONAL && lp.terms()>11) {
			if(PARSE_OPTIONAL_MD_ONLY){
				optional=new ArrayList<String>(1);
				for(int i=11, terms=lp.terms(); i<terms; i++) {
					if(lp.termStartsWith("MD:", i)){
						String s=lp.parseString(i);
						optional.add(s);
						lp.incrementA(3);
						mdTag=lp.parseByteArrayFromCurrentField();
					}
				}
			}else{
				optional=new ArrayList<String>(lp.terms()-11);
				for(int i=11, terms=lp.terms(); i<terms; i++) {
					String s=lp.parseString(i);
					optional.add(s);
				}
			}
		}
		
		trimNames();
	}
	
	//TODO: Should be replaced by LineParser
	public SamLine(byte[] s){
		assert(s[0]!='@') : "Tried to make a SamLine from a header: "+new String(s);
		
		int a=0, b=0;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(s);
		if(PARSE_0){qname=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(s);
		flag=Parse.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(s);
		if(PARSE_2){
			if(RNAME_AS_BYTES){
				rname=(b==a+1 && s[a]=='*' ? null : KillSwitch.copyOfRange(s, a, b));
			}else{
				rnameS=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));
			}
		}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(s);
		pos=Parse.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(s);
		mapq=Parse.parseInt(s, a, b);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(s);
		if(PARSE_5){cigar=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(s);
		if(PARSE_6){rnext=(b==a+1 && s[a]=='*' ? null : KillSwitch.copyOfRange(s, a, b));}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(s);
		if(PARSE_7){pnext=(b==a+1 && s[a]=='*' ? 0 :Parse.parseInt(s, a, b));}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(s);
		if(PARSE_8){tlen=Parse.parseInt(s, a, b);}
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 9: "+new String(s);
//		seq=new String(s, a, b-a);
		seq=(b==a+1 && s[a]=='*' ? null : KillSwitch.copyOfRange(s, a, b));
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 10: "+new String(s);
//		qual=new String(s, a, b-a);
		if(PARSE_10){qual=(b==a+1 && s[a]=='*' ? null : KillSwitch.copyOfRange(s, a, b));}
		b++;
		a=b;

		assert((seq==bytestar)==(Tools.equals(seq, bytestar)));
		assert((qual==bytestar)==(Tools.equals(qual, bytestar)));
		
		if(mapped() && strand()==Shared.MINUS && FLIP_ON_LOAD){
			if(seq!=bytestar){AminoAcid.reverseComplementBasesInPlace(seq);}
			if(qual!=bytestar){Tools.reverseInPlace(qual);}
		}
		
		if(qual!=null && qual!=bytestar){
			for(int i=0; i<qual.length; i++){qual[i]-=33;}
		}

		if(PARSE_OPTIONAL && b<s.length) {
			if(PARSE_OPTIONAL_MD_ONLY) {
				optional=new ArrayList<String>(1);
				while(b<s.length){
					while(b<s.length && s[b]!='\t'){b++;}
					if(b>a){
						if(b>=a+3 && s[a]=='M' && s[a+1]=='D' && s[a+2]==':'){
							String x=new String(s, a, b-a);
							optional.add(x);
							return;
						}
					}else{
						//Empty field
					}
					b++;
					a=b;
				}
			}else if(PARSE_OPTIONAL_MATEQ_ONLY) {
				optional=new ArrayList<String>(1);
				while(b<s.length){
					while(b<s.length && s[b]!='\t'){b++;}
					if(b>a){
						if(b>=a+3 && s[a]=='Y' && s[a+1]=='Q' && s[a+2]==':'){
							String x=new String(s, a, b-a);
							optional.add(x);
							return;
						}
					}else{
						//Empty field
					}
					b++;
					a=b;
				}
			}else{
				optional=new ArrayList<String>(4);
				while(b<s.length){
					while(b<s.length && s[b]!='\t'){b++;}
					if(b>a){
						String x=new String(s, a, b-a);
						optional.add(x);
					}else{
						//Empty field
					}
					b++;
					a=b;
				}
			}
		}
		
		trimNames();
	}
	
	public void trimNames(){
//		System.err.println();
//		System.err.println("rname= "+new String(rname));
//		System.err.println("qname= "+new String(qname));
		if(Shared.TRIM_RNAME){
			if(RNAME_AS_BYTES){
				setRname(Tools.trimToWhitespace(rname()));
				setRnext(Tools.trimToWhitespace(rnext()));
			}else{
				setRname(Tools.trimToWhitespace(rnameS()));
				setRnext(Tools.trimToWhitespace(rnext()));
			}
		}
		if(Shared.TRIM_READ_COMMENTS){
			qname=(Tools.trimToWhitespace(qname));
		}
//		System.err.println("rname2="+new String(rname));
//		System.err.println("qname2="+new String(qname));
//		assert(false) : Shared.TRIM_RNAME+", "+Shared.TRIM_READ_COMMENTS+", "+new String(rname)+", "+qname;
	}
	
	public static final int parseFlagOnly(byte[] s){
		assert(s!=null && s.length>0) : "Blank line.";
		if(s[0]=='@'){return -1;}
		
		int a=0, b=0;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(s);
		b++;
		a=b;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(s);
		int flag=Parse.parseInt(s, a, b);
		return flag;
	}
	
	public static final String parseNameOnly(byte[] s){
		assert(s!=null && s.length>0) : "Blank line.";
		if(s[0]=='@'){return null;}
		
		int a=0, b=0;
		
		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(s);
		String qname=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a));
		return qname;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Cigar            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static String toCigar13(byte[] match, int readStart, int readStop, long reflen, byte[] bases){
		if(match==null || readStart==readStop){return null;}
		ByteBuilder sb=new ByteBuilder(8);
		int count=0;
		char mode='=';
		char lastMode='=';
		
		int refloc=readStart;

		int cigarlen=0; //for debugging
		int opcount=0; //for debugging
		
		for(int mpos=0; mpos<match.length; mpos++){

			byte m=match[mpos];
			
			boolean sfdflag=false;
			if(SOFT_CLIP && (refloc<0 || refloc>=reflen)){
				mode='S'; //soft-clip out-of-bounds
				if(m!='I'){refloc++;}
				if(m=='D'){sfdflag=true;} //Don't add soft-clip count for deletions!
			}else if(m=='m' || m=='s' || m=='S' || m=='N' || m=='B'){//Little 's' is for a match classified as a sub to improve the affine score.
				mode='M';
				refloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				refloc++;
			}else if(m=='C'){
				mode='S';
				refloc++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}

			if(mode!=lastMode){
				if(count>0){//Prevents an initial length-0 match
					sb.append(count);
//					sb.append(lastMode);
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}
					else{sb.append(lastMode);}
					if(lastMode!='D'){cigarlen+=count;}
					opcount+=count;
				}
				count=0;
				lastMode=mode;
			}

			count++;
			if(sfdflag){count--;}
		}
		sb.append(count);
		if(mode=='D' && count>INTRON_LIMIT){sb.append('N');}
		else{sb.append(mode);}
		if(mode!='D'){cigarlen+=count;}
		opcount+=count;
		
		assert(bases==null || cigarlen==bases.length) : "\n(cigarlen = "+cigarlen+") != (bases.length = "+(bases==null ? -1 : bases.length)+")\n" +
				"cigar = "+sb+"\nmatch = "+new String(match)+"\nbases = "+new String(bases)+"\n";
		
		return sb.toString();
	}
	
	public static String toCigar13(String cigar14) {
		if(cigar14==null){return null;}
		final int len=cigar14.length();

		int current=0;
		int mcount=0;
		ByteBuilder sb=new ByteBuilder(len);
		
		for(int i=0; i<len; i++){
			char b=cigar14.charAt(i);
			if(Tools.isDigit(b)){
				current=(10*current)+(b-'0');
			}else{
				if(b=='X' || b=='=' || b=='M'){
					mcount+=current;
				}else{
					if(mcount>0){
						sb.append(mcount);
						sb.append('M');
						mcount=0;
					}
					sb.append(current);
					sb.append(b);
				}
				current=0;
			}
		}
		assert(current==0);
		if(mcount>0){
			sb.append(mcount);
			sb.append('M');
			mcount=0;
		}
		return sb.toString();
	}
	
	
	public static String toCigar14(byte[] match, int readStart, int readStop, long reflen, byte[] bases){
//		assert(false) : readStart+", "+readStop+", "+reflen;
		if(match==null || readStart==readStop){return null;}
		ByteBuilder sb=new ByteBuilder(8);
		int count=0;
		char mode='=';
		char lastMode='=';
		
		int refloc=readStart;

		int cigarlen=0; //for debugging
		int opcount=0; //for debugging
		
		for(int mpos=0; mpos<match.length; mpos++){

			byte m=match[mpos];
			
			boolean sfdflag=false;
			if(SOFT_CLIP && (refloc<0 || refloc>=reflen)){
				mode='S'; //soft-clip out-of-bounds
				if(m!='I'){refloc++;}
				if(m=='D'){sfdflag=true;} //Don't add soft-clip count for deletions!
			}else if(m=='m' || m=='s'){//Little 's' is for a match classified as a sub to improve the affine score.
				mode='=';
				refloc++;
			}else if(m=='S' || m=='V'){
				mode='X';
				refloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				refloc++;
			}else if(m=='C'){
				mode='S';
				refloc++;
			}else if(m=='N' || m=='B'){
				mode='M';
				refloc++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}

			if(mode!=lastMode){
				if(count>0){//Prevents an initial length-0 match
					sb.append(count);
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}
					else{sb.append(lastMode);}
					if(lastMode!='D'){cigarlen+=count;}
					opcount+=count;
				}
				count=0;
				lastMode=mode;
			}

			count++;
			if(sfdflag){count--;}
		}
		sb.append(count);
		if(mode=='D' && count>INTRON_LIMIT){
			sb.append('N');
		}else{
			sb.append(mode);
		}
		if(mode!='D'){cigarlen+=count;}
		opcount+=count;
		
		assert(bases==null || cigarlen==bases.length) : "\n(cigarlen = "+cigarlen+") != (bases.length = "+(bases==null ? -1 : bases.length)+")\n" +
				"cigar = "+sb+"\nmatch = "+new String(match)+"\nbases = "+new String(bases)+"\n";
		
		return sb.toString();
	}

	public boolean cigarContainsOnlyME() {
		if(cigar==null || cigar.length()==0){return false;}
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c) || c=='M' || c=='='){
				//do nothing;
			}else{
				return false;
			}
		}
		return true;
	}
	
	public int calcCigarLength(boolean includeSoftClip, boolean includeHardClip){
		return calcCigarLength(cigar, includeSoftClip, includeHardClip);
	}
	
	public int calcCigarReadLength(boolean includeSoftClip, boolean includeHardClip){
		return calcCigarReadLength(cigar, includeSoftClip, includeHardClip);
	}
	
	/** Reference length of cigar string */
	public static int calcCigarLength(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='=' || c=='X' || c=='D' || c=='N'){
					len+=current;
				}else if(c=='S'){
					if(includeSoftClip){len+=current;}
				}else if (c=='H'){
					//In this case, the base string is the wrong length since letters were truncated.
					//Therefore, the bases cannot be used for calling variations after mapping.
					//Hard clipping messes up original location verification.
					//Therefore...  len+=current would be best in practice, but for GRADING purposes, leaving it disabled is best.

					if(includeHardClip){len+=current;}
				}else if(c=='I'){
					//do nothing
				}else if(c=='P'){
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
					//'P' is currently poorly defined
				}else{
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
				}
				current=0;
			}
		}
		return len;
	}
	
	/** Reference length of cigar string */
	public static int calcCigarReadLength(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='=' || c=='X' || c=='I'){
					len+=current;
				}else if(c=='S'){
					if(includeSoftClip){len+=current;}
				}else if (c=='H'){
					//In this case, the base string is the wrong length since letters were truncated.
					//Therefore, the bases cannot be used for calling variations after mapping.
					//Hard clipping messes up original location verification.
					//Therefore...  len+=current would be best in practice, but for GRADING purposes, leaving it disabled is best.

					if(includeHardClip){len+=current;}
				}else if(c=='D' || c=='N'){
					//do nothing
				}else if(c=='P'){
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
					//'P' is currently poorly defined
				}else{
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
				}
				current=0;
			}
		}
		return len;
	}
	
	/** Number of query bases in cigar string */
	public static int calcCigarBases(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='=' || c=='X' || c=='I'){
					len+=current;
				}else if(c=='D' || c=='N'){
					//do nothing
				}else if (c=='H'){
					if(includeHardClip){len+=current;}
				}else if(c=='S'){
					if(includeSoftClip){len+=current;}
				}else if(c=='P'){
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
					//'P' is currently poorly defined
				}else{
					throw new RuntimeException("Unhandled cigar symbol: "+c+"\n"+cigar+"\n");
				}
				current=0;
			}
		}
		return len;
	}
	
	/** Length of clipped initial bases.  Used to calculate correct start location of clipped reads. */
	public static int countLeadingClip(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null || (!includeSoftClip && !includeHardClip)){return 0;}
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isLetter(c) || c=='='){
				if(c=='H'){
					if(includeHardClip){
						len+=current;
					}
				}else if(c=='S'){
					if(includeSoftClip){
						len+=current;
					}
				}else{
					break;
				}
				current=0;
			}else{
				current=(current*10)+(c-'0');
			}
		}
		return len;
	}
	
	/** Length of clipped final bases.  Used to calculate correct stop location of clipped reads. */
	public static int countTrailingClip(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null || (!includeSoftClip && !includeHardClip)){return 0;}
		int len=0;
		if(includeHardClip){len+=countTrailingHardClip(cigar);}
		int last=cigar.lastIndexOf('S');
		
		int mult=1;
		int i;
		for(i=last-1; i>=0; i--){
			char c=cigar.charAt(i);
			if(Tools.isLetter(c) || c=='='){
				break;
			}
			len+=(len+(c-'0')*mult);
			mult*=10;
		}
		if(i<0){return 0;}
		return len;
	}
	
	/** Length of clipped final bases.  Used to calculate correct stop location of clipped reads. */
	public static int countTrailingHardClip(String cigar){
		if(cigar==null){return 0;}
		int last=cigar.lastIndexOf('H');
		
		int mult=1, len=0;
		int i;
		for(i=last-1; i>=0; i--){
			char c=cigar.charAt(i);
			if(Tools.isLetter(c) || c=='='){
				break;
			}
			len+=(len+(c-'0')*mult);
			mult*=10;
		}
		if(i<0){return 0;}
		return len;
	}
	
	public static int countMdSubs(String mdTag){
		assert(mdTag!=null);

		final int NORMAL=0, SUB=1, DEL=2;
		int dels=0, subs=0, normals=0;
		
		if(mdTag!=null){
			int current=0;
			int mode=NORMAL;
			int i=0;
			if(mdTag.startsWith("MD:Z:")){i=5;}
			for(final int max=mdTag.length(); i<max; i++){
				char c=mdTag.charAt(i);
				if(Tools.isDigit(c)){
					current=(current*10)+(c-'0');
					mode=NORMAL;
				}else{
					if(current>0){
						if(mode==NORMAL){normals+=current;}
						else{assert(false) : mode+", "+current;}
						current=0;
					}
					if(c=='^'){mode=DEL;}
					else if(mode==DEL){
						dels++;
					}else if(mode==NORMAL || mode==SUB){
						mode=SUB;
						subs++;
					}
				}
			}
		}
		return subs;
	}
	
	/** Length of clipped initial bases. */
	public static int countLeadingClip(byte[] match){
		if(match==null || match.length<1 || match[0]!='C'){return 0;}
		int clips=0;
		int current=0;
		for(int mloc=0; mloc<match.length; mloc++){
			byte b=match[mloc];
			if(Tools.isDigit(b)){
				current=current*10+(b-'0');
			}else{
				if(current>0){
					clips=clips+current-1;
				}
				current=0;
				if(b!='C'){break;}
				clips++;
			}
		}
		if(current>0){
			clips=clips+current-1;
		}
		return clips;
	}
	
	/** Length of match string portion describing clipped initial bases. */
	public static int countLeadingClip2(byte[] match){
		if(match==null || match.length<1 || match[0]!='C'){return 0;}
		int mloc=0;
		for(; mloc<match.length; mloc++){
			byte b=match[mloc];
			if(b!='C' && !Tools.isDigit(b)){return mloc;}
		}
		return match.length;
	}
	
	/** Length of clipped trailing bases. */
	public static int countTrailingClip(byte[] match){
		if(match==null){return 0;}
		int clips=0;
		for(int mloc=match.length-1; mloc>=0; mloc--){
			byte b=match[mloc];
			assert(!Tools.isDigit(b)) : new String(match);
			if(b=='C'){
				clips++;
			}else{
				break;
			}
		}
		return clips;
	}
	
	/** Length of clipped (out of bounds) initial insertions and deletions. */
	public static int countLeadingIndels(int rloc, byte[] match){
		if(match==null || rloc>=0){return 0;}
		int dels=0;
		int inss=0;
		int cloc=0;
		for(int mloc=0; mloc<match.length && rloc<0; mloc++){
			byte b=match[mloc];
			assert(!Tools.isDigit(b));
			if(b=='D'){
				dels++;
				rloc++;
			}else if(b=='I'){
				inss++;
				cloc++;
			}else{
				rloc++;
				cloc++;
			}
		}
		return dels-inss;
	}
	
	/** Length of clipped (out of bounds) trialing insertions and deletions. */
	public static int countTrailingIndels(int rloc, int rlen, byte[] match){
		if(match==null || rloc>=0){return 0;}
		int dels=0;
		int inss=0;
		int cloc=0;
		for(int mloc=match.length; mloc>=0 && rloc>=rlen; mloc--){
			byte b=match[mloc];
			assert(!Tools.isDigit(b));
			if(b=='D'){
				dels++;
				rloc--;
			}else if(b=='I'){
				inss++;
				cloc--;
			}else{
				rloc--;
				cloc--;
			}
		}
		return dels-inss;
	}

	public int mappedNonClippedBases() {
		if(!mapped() || cigar==null){return 0;}
		
		int len=0;
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='='){
					len+=current;
				}else if(c=='X'){
					len+=current;
				}else if(c=='D' || c=='N'){
					
				}else if(c=='I'){
					len+=current;
				}else if(c=='S' || c=='H' || c=='P'){
					
				}
				current=0;
			}
		}
		return len;
	}
	
	/**
	 * @param cigar
	 * @return Max consecutive match, sub, del, ins, or clip symbols
	 */
	public static final int[] cigarToMdsiMax(String cigar) {
		if(cigar==null){return null;}
		int[] msdic=KillSwitch.allocInt1D(5);
		
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='='){
					msdic[0]=Tools.max(msdic[0], current);
				}else if(c=='X'){
					msdic[1]=Tools.max(msdic[1], current);
				}else if(c=='D' || c=='N'){
					msdic[2]=Tools.max(msdic[2], current);
				}else if(c=='I'){
					msdic[3]=Tools.max(msdic[3], current);
				}else if(c=='S' || c=='H' || c=='P'){
					msdic[4]=Tools.max(msdic[4], current);
				}
				current=0;
			}
		}
		return msdic;
	}

	public float calcIdentity() {
		assert(cigar!=null);
		int match=0, other=0;
		
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='='){
					match+=current;
				}else if(c=='M'){
					
				}else if(c=='X'){
					other+=current;
				}else if(c=='D'){
					other+=current;
				}else if(c=='N'){
					
				}else if(c=='I'){
					other+=current;
				}else if(c=='S' || c=='H' || c=='P'){
					
				}
				current=0;
			}
		}
		return match/(float)Tools.max(match+other, 1);
	}
	
	public int countSubs() {
		if(cigar==null){return 0;}
		
		int current=0;
		int subs=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='X'){
					subs+=current;
				}
				current=0;
			}
		}
		return subs;
	}
	
	/**
	 * @param cigar
	 * @return Total number of match, sub, del, ins, or clip symbols
	 */
	public static final int[] cigarToMsdic(String cigar) {
		if(cigar==null){return null;}
		int[] msdic=KillSwitch.allocInt1D(5);
		
		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='M' || c=='='){
					msdic[0]+=current;
				}else if(c=='X'){
					msdic[1]+=current;
				}else if(c=='D' || c=='N'){
					msdic[2]+=current;
				}else if(c=='I'){
					msdic[3]+=current;
				}else if(c=='S' || c=='H' || c=='P'){
					msdic[4]+=current;
				}
				current=0;
			}
		}
		return msdic;
	}
	
	/**
	 * @param allowM Allow M symbols in the cigar string
	 * @return Match string of this cigar string when possible, otherwise null.
	 * Takes into account MD tag and bases, but not reference (other than in MD tag).
	 */
	public final byte[] toShortMatch(boolean allowM) {
		if(cigar==null || cigar.equals(stringstar)){return null;}
		if(allowM){return cigarToShortMatch_old(cigar, allowM);}
		
//		System.err.println("\nInput: cigar="+cigar+", MD="+mdTag());//123

		final boolean fixMatchSubs;
		final boolean fixMatchNs;
		boolean foundE=false;
		boolean foundX=false;
		boolean foundM=false;
//		System.err.println("Block 1.");//123
		{
			int current=0, total=0;
			for(int i=0; i<cigar.length(); i++){
				char c=cigar.charAt(i);

				if(Tools.isDigit(c)){
					current=(current*10)+(c-'0');
				}else{

					if(c=='H'){
						current=0; //Information destroyed
					}else if(c=='P'){
						return null; //Undefined symbol
					}
					foundE|=(c=='=');
					foundX|=(c=='X');
					foundM|=(c=='M');

					total+=current;
					current=0;
				}
			}
			if(total<1){return null;}
			fixMatchSubs=(!allowM && foundM && !foundX && !foundE);//Note: allowM already exited.
			fixMatchNs=(FIX_MATCH_NS && foundX && !foundM); //Means no-calls are possibly marked as X, which is technically OK.

//			System.err.println("allowM="+allowM);//123
//			System.err.println("foundE="+foundE);//123
//			System.err.println("foundX="+foundX);//123
//			System.err.println("foundM="+foundM);//123
		}
		
//		System.err.println("Block 2.");//123
		final String mdTag;
		final int mdSubs;
		final byte[] refBases;
		
		//1) if fixMatch, grab MD tag
		//2) if MD and no subs, return
		//3) grab ref bases
		
		if(fixMatchSubs || fixMatchNs){
			final String md0=mdTag();
			mdSubs=(md0==null ? -1 : countMdSubs(md0));
			if(mdSubs==0 && !fixMatchNs){
				refBases=null;
				mdTag=null;
			}else{
				mdTag=md0;
				if(mdTag!=null && PREFER_MDTAG){
					refBases=null;
				}else{
					ScafMap map=ScafMap.defaultScafMap();
					assert(!fixMatchSubs || mdTag!=null || map!=null) : "TODO: Encountered a read with 'M' in cigar string but no MD tag and no ScafMap loaded.\n"
							+ "This can normally be resolved by adding the flag ref=file, where file is the fasta file to which the reads were mapped.\n\n"+this;
					if(map==null){
						refBases=null;
					}else{
						Scaffold scaf=map.getScaffold(rnameS());
						assert(!fixMatchSubs || mdTag!=null || scaf!=null) : "Encountered a read with 'M' in cigar string but no scaffold loaded for "+rnameS();
						if(scaf==null){
							refBases=null;
						}else{
							refBases=scaf.bases;
							assert(!fixMatchSubs || mdTag!=null || refBases!=null) : "Encountered a read with 'M' in cigar string but no ref bases loaded for "+rnameS();
							if(fixMatchSubs && refBases==null && mdTag==null){
								return null;
							}
						}
					}
				}
			}
		}else{
			mdSubs=-1;
			mdTag=null;
			refBases=null;
		}
		
//		System.err.println("mdTag="+mdTag);//123
//		System.err.println("mdSubs="+mdSubs);//123
		
		char mSymbol=((foundX || foundE || !foundM) ? 'N' : mdSubs>=0 ? 'm' : 'N');
		
//		System.err.println("Block 3.");//123
		final byte[] match0;
		{
			ByteBuilder sb=new ByteBuilder(cigar.length());
			int current=0;
			for(int cpos=0, max=cigar.length(); cpos<max; cpos++){
				char c=cigar.charAt(cpos);
				if(Tools.isDigit(c)){
					current=(current*10)+(c-'0');
				}else{
					if(c=='='){
						sb.append('m');
						if(current>1){sb.append(current);}
					}else if(c=='X'){
						sb.append('S');
						if(current>1){sb.append(current);}
					}else if(c=='D' || c=='N'){
						sb.append('D');
						if(current>1){sb.append(current);}
					}else if(c=='I'){
						sb.append('I');
						if(current>1){sb.append(current);}
					}else if(c=='S'){
						sb.append('C');
						if(current>1){sb.append(current);}
					}else if(c=='M'){
						sb.append(mSymbol);
						if(current>1){sb.append(current);}
					}
					current=0;
				}
			}
			match0=(sb.array.length==sb.length() ? sb.array : sb.toBytes());
//			System.err.println("match="+new String(match0));//123
		}
		
//		System.err.println("Block 4.");//123
		
		if((!fixMatchSubs || mdSubs==0) && (!fixMatchNs || seq==null || refBases==null)){return match0;}
		assert(((mdTag!=null || refBases!=null) && fixMatchSubs && mdSubs!=0) || fixMatchNs /*|| noCalls>0*/) : 
			(mdTag!=null)+", "+(refBases!=null)+", "+(fixMatchSubs)+", "+(mdSubs!=0)+", "+fixMatchNs/*+", "+(noCalls)*/;
		
//		assert(false) : mdTag+", "+refBases+", "+processMD+", "+mdSubs+"\n"+this;
		
//		System.err.println("processMD="+processMD+", mdSubs="+mdSubs+", mdTag="+mdTag);//123
		
		final byte[] bases;
		if(refBases!=null){
//			bases=(strand()==1) ? AminoAcid.reverseComplementBases(seq) : seq;//Why not reverse in place?
			if(strand()==1){AminoAcid.reverseComplementBasesInPlace(seq);}
			bases=seq;
		}else{bases=null;}

//		System.err.println("Block 5.");//123
		final byte[] longmatch=Read.toLongMatchString(match0);

//		final int noCalls=((foundE || foundX) && foundM ? 1 : seq==null ? -1 : Read.countNocalls(seq));
		
//		System.err.println("Block 6");//123
		
		if(mdTag!=null && (refBases==null || PREFER_MDTAG)){
//			System.err.println("match="+new String(longmatch));//123

			final int noCalls=seq==null ? -1 : Read.countNocalls(seq);
			if(noCalls>0 && bases!=null && refBases==null){//Not necessary if ref sequence is present
				int bpos=0;
				for(int mpos=0; mpos<longmatch.length; mpos++){
					final byte m=longmatch[mpos];
					if(m=='C'){
						bpos++;
					}else if(m=='m' || m=='s' || m=='S' || m=='N'){
						if(!AminoAcid.isFullyDefined(bases[bpos])){longmatch[mpos]='N';}
						bpos++;
					}else if(m=='I' || m=='X' || m=='Y'){
						bpos++;
					}else if(m=='D'){
						//do nothing
					}else{
						assert(false) : m;
					}
				}
			}
			
			final MDWalker walker=new MDWalker(mdTag, cigar, longmatch, this);
			
			walker.fixMatch(bases);
		}else
		if(refBases!=null){
			final int refStart=start(true, false);
			fixMatch(bases, refBases, longmatch, refStart, false);
		}
		
//		else if(refBases!=null){
//			final int refStart=start(true, false);
//			fixMatch(bases, refBases, longmatch, refStart, false);
//		}
		else{
			assert(false) : "Fallthorugh.";
		}
		
		final byte[] match=Read.toShortMatchString(longmatch);
		
		if(bases!=null && strand()==1){AminoAcid.reverseComplementBasesInPlace(seq);}
		
//		System.err.println("Block 7.");//123
//		System.err.println("Returning "+new String(match));//123
		return match;
	}
	
	/** Requires longmatch.
	 * Replaces M  */
	public static void fixMatch(byte[] call, byte[] ref, byte[] match, int refstart, boolean unClip){
		for(int mpos=0, rpos=refstart, cpos=0; mpos<match.length; mpos++){
			assert(cpos>=0 && cpos<call.length) : "\n"+new String(match)+"\n"+new String(call)+"\n"+mpos+", "+cpos;
			final byte m=match[mpos];
			
			if(rpos<0 || rpos>=ref.length){
				if(m=='I'){
					assert(false) : "Insertion off scaffold end: "+refstart+", "+ref.length+"\n"+new String(call)+"\n"+new String(match);
					cpos++;
				}else if(m=='D'){
					assert(false) : "Deletion off scaffold end: "+refstart+", "+ref.length+"\n"+new String(call)+"\n"+new String(match);
					rpos++;	
				}else{
					match[mpos]='C';
					rpos++;
					cpos++;
				}
			}else if(m=='m' || m=='S' || m=='N' || m=='s' || (m=='C' && unClip)){
				final byte c=Tools.toUpperCase(call[cpos]);
				final byte r=Tools.toUpperCase(ref[rpos]);
				final boolean defined=(AminoAcid.isFullyDefined(c) && AminoAcid.isFullyDefined(r));
				if(!defined){
					match[mpos]='N';
				}else if(c==r){
					match[mpos]='m';
				}else{
					match[mpos]='S';
				}
				rpos++;
				cpos++;
			}else if(m=='C'){ //Do nothing for clipped call
				rpos++;
				cpos++;
			}else if(m=='I' || m=='X' || m=='Y'){
				cpos++;
			}else if(m=='D'){
				rpos++;
			}else{
				assert(false) : Character.toString((char)m);
			}
		}
	}
	
	/**
	 * @param cigar
	 * @return Match string of this cigar string when possible, otherwise null
	 */
	public static final byte[] cigarToShortMatch_old(String cigar, boolean allowM) {
		
		int current=0;
		ByteBuilder sb=new ByteBuilder(cigar.length());
		
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='='){
					sb.append('m');
					if(current>1){sb.append(current);}
				}else if(c=='X'){
					sb.append('S');
					if(current>1){sb.append(current);}
				}else if(c=='D' || c=='N'){
					sb.append('D');
					if(current>1){sb.append(current);}
				}else if(c=='I'){
					sb.append('I');
					if(current>1){sb.append(current);}
				}else if(c=='S'){
					sb.append('C');
					if(current>1){sb.append(current);}
				}else if(c=='M'){
					if(!allowM){return null;}
//					sb.append('B');
					sb.append('N');
					if(current>1){sb.append(current);}
				}
				current=0;
			}
		}
		
		if(sb.array.length==sb.length()){return sb.array;}
		return sb.toBytes();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Tags             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static String makeStopTag(int pos, int seqLength, String cigar, boolean perfect){
//		return "YS:i:"+(pos+((cigar==null || perfect) ? seqLength : -countLeadingClip(cigar, false)+calcCigarLength(cigar, false))-1); //123456789
		return "YS:i:"+(pos+((cigar==null || perfect) ? seqLength : calcCigarLength(cigar, true, false))-1);
	}
	
	public static String makeLengthTag(int pos, int seqLength, String cigar, boolean perfect){
		if(cigar==null || perfect){return "YL:Z:"+seqLength+","+seqLength;}
		return "YL:Z:"+(seqLength-countLeadingClip(cigar, true, false))+","+calcCigarLength(cigar, false, false);
	}
	
	public static String makeIdentityTag(byte[] match, boolean perfect){
		if(perfect){return "YI:f:100";}
		float f=Read.identity(match);
		return Tools.format("YI:f:%.2f", (100*f));
	}
	
	public static String makeScoreTag(int score){
		return "YR:i:"+score;
	}
	
	public String matchTag(){
		if(optional==null){return null;}
		for(String s : optional){
			if(s.startsWith("X2:Z:")){
				return s;
			}
		}
		return null;
	}
	
	private String makeXSTag(Read r){
		if(r.mapped() && cigar!=null && cigar.indexOf('N')>=0){
//			System.err.println("For read "+r.pairnum()+" mapped to strand "+r.strand());
			boolean plus=(r.strand()==Shared.PLUS); //Assumes secondstrand=false
//			System.err.println("plus="+plus);
			if(r.pairnum()!=0){plus=!plus;}
//			System.err.println("plus="+plus);
			if(XS_SECONDSTRAND){plus=!plus;}
//			System.err.println("plus="+plus);
			return (plus ? XSPLUS : XSMINUS);
		}else{
			return null;
		}
	}
	
	public static String makeMdTag(int chrom, int refstart, byte[] match, byte[] call, int scafloc, int scaflen){
		if(match==null || chrom<0){return null;}
		ByteBuilder md=new ByteBuilder(8);
		md.append("MD:Z:");
		
		ChromosomeArray cha=Data.getChromosome(chrom);
		
		final int scafstop=scafloc+scaflen;
		
		byte prevM='?';
		int count=0;
		int dels=0;
		boolean prevSub=false;
		for(int mpos=0, rpos=refstart, cpos=0; mpos<match.length; mpos++){
			assert(cpos>=0 && cpos<call.length) : "\n"+new String(match)+"\n"+new String(call)+"\n"+mpos+", "+cpos+", "+dels+", "+INTRON_LIMIT;
			final byte c=call[cpos];
			final byte m=match[mpos];
			
			if(prevM=='D' && m!='D'){
				if(dels<=INTRON_LIMIT){//Otherwise, ignore it
					md.append(count);
					count=0;
					md.append('^');
					for(int i=rpos-dels; i<rpos; i++){
						md.append((char)cha.get(i));
					}
					dels=0;
				}
			}
			
			if(m=='C' || rpos<scafloc || rpos>=scafstop){ //Do nothing for clipped bases
				rpos++;
				if(m!='D'){cpos++;}
			}else if(m=='m' || m=='s'){
				count++;
				rpos++;
				cpos++;
			}else if(m=='S'){
				if(count>0 || !prevSub){md.append(count);}
				md.append((char)cha.get(rpos));

				count=0;
				rpos++;
				cpos++;
				prevSub=true;
			}else if(m=='N'){
				
				final byte r=cha.get(rpos);
				
				if(c==r){//Act like match
					count++;
					rpos++;
					cpos++;
				}else{//Act like sub
					if(count>0 || !prevSub){md.append(count);}
					md.append((char)r);

					count=0;
					rpos++;
					cpos++;
					prevSub=true;
				}
			}else if(m=='I' || m=='X' || m=='Y'){
				cpos++;
//				count++;
			}else if(m=='D'){
//				if(prevM!='D'){
//					md.append(count);
//					count=0;
//					md.append('^');
//				}
//				md.append((char)cha.get(rpos));
				
				rpos++;
				dels++;
			}
			prevM=m;
			
		}
//		if(count>0){
			md.append(count);
//		}
		
		return md.toString();
	}
	
	public static int calcLeftClip(String cig, String id){
		if(cig==null){return 0;}
		int len=0;
		for(int i=0; i<cig.length(); i++){
			char c=cig.charAt(i);
			if(Tools.isDigit(c)){
				len=len*10+(c-'0');
			}else{
				assert(c!='S' || i<cig.length()-1);//ban entirely soft-clipped reads
				return (c=='S') ? len : 0;
			}
		}
		return 0;
	}
	
	public static int calcRightClip(String cig, String id){
		if(cig==null || cig.length()<1 || cig.charAt(cig.length()-1)!='S'){return 0;}
		int pos=cig.length()-2;
		for(; pos>=0 && Tools.isDigit(cig.charAt(pos)); pos--){}
		
		assert(pos>0) : cig+", id="+id+", pos="+pos;//ban entirely soft-clipped reads
		
		int len=0;
		for(int i=pos+1; i<cig.length(); i++){
			char c=cig.charAt(i);
			if(Tools.isDigit(c)){
				len=len*10+(c-'0');
			}else{
				return (c=='S') ? len : 0;
			}
		}
		return len;
	}
	
	public ArrayList<String> makeOptionalTags(Read r, Read r2, boolean perfect, int scafloc, int scaflen, boolean inbounds, boolean inbounds2){
		if(NO_TAGS){return null;}
		final boolean mapped=r.mapped();
		if(!mapped && READGROUP_ID==null && !MAKE_CUSTOM_TAGS && !MAKE_TIME_TAG){return null;}

		ArrayList<String> optionalTags=new ArrayList<String>(8);

		if(mapped){
			if(!r.secondary() && r.ambiguous() && MAKE_XT_TAG){optionalTags.add("XT:A:R");} //Not sure what do do for secondary alignments

//			int nm=r.length();
//			int dels=0;
			
			int nm=0;
			
//			//Only works for cigar strings in format 1.4+
//			if(perfect){nm=0;}else if(cigar!=null){
//				int len=0;
//				for(int i=0; i<cigar.length(); i++){
//					char c=cigar.charAt(i);
//					if(Tools.isDigit(c)){
//						len=len*10+(c-'0');
//					}else{
//						if(c=='X' || c=='I' || c=='D' || c=='M'){
//							nm+=len;
//						}
//						len=0;
//					}
//				}
////				System.err.println("\nRead "+r.id+": nm="+nm+"\n"+cigar+"\n"+new String(r.match));
//				System.err.println("\nRead "+r.id+": nm="+nm);
//			}
			
			if(perfect){nm=0;}else if(r.match!=null){
				nm=0;
				int leftclip=calcLeftClip(cigar, r.id), rightclip=calcRightClip(cigar, r.id);
				final int from=leftclip, to=r.length()-rightclip;
				int delsCurrent=0;
				for(int i=0, cpos=0; i<r.match.length; i++){
					final byte b=r.match[i];
					
//					System.err.println("i="+i+", cpos="+cpos+", from="+from+", ")
					
					if(cpos>=from && cpos<to){
						if(b=='I' || b=='S' || b=='N' || b=='X' || b=='Y'){nm++;}

						if(b=='D'){delsCurrent++;}
						else{
							if(delsCurrent<=INTRON_LIMIT){nm+=delsCurrent;}
							delsCurrent=0;
						}
					}
					if(b!='D'){cpos++;}
				}
				if(delsCurrent<=INTRON_LIMIT){nm+=delsCurrent;}
				//				assert(false) : nm+",  "+dels+", "+delsCurrent+", "+r.length()+", "+r.match.length;

//				assert(false) : "rlen="+r.length()+", nm="+nm+", dels="+delsCurrent+", intron="+INTRON_LIMIT+", inbound1="+inbounds+", ib2="+inbounds2+"\n"+new String(r.match);
				
//				System.err.println("\nRead "+r.id+": left="+leftclip+", right="+rightclip+", nm="+nm+"\n"+cigar+"\n"+new String(r.match));
				
			}
			
			if(MAKE_NM_TAG){
				if(perfect){optionalTags.add("NM:i:0");}
				else if(r.match!=null){optionalTags.add("NM:i:"+(nm));}
			}
			if(MAKE_SM_TAG){optionalTags.add("SM:i:"+mapq);}
			if(MAKE_AM_TAG){optionalTags.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, r2.mapScore/r2.length()) : 0)));}
			
			if(MAKE_TOPHAT_TAGS){
				optionalTags.add("AS:i:0");
				if(cigar==null || cigar.indexOf('N')<0){
					optionalTags.add("XN:i:0");
				}else{
				}
				optionalTags.add("XM:i:0");
				optionalTags.add("XO:i:0");
				optionalTags.add("XG:i:0");
				if(cigar==null || cigar.indexOf('N')<0){
					optionalTags.add("YT:Z:UU");
				}else{
				}
				optionalTags.add("NH:i:1");
			}else if(MAKE_XM_TAG){//XM tag.  For bowtie compatibility; unfortunately it is poorly defined.
				int x=0;
				if(r.discarded() || (!r.ambiguous() && !mapped)){
					x=0;//TODO: See if the flag needs to be present in this case.
				}else if(mapped){
					x=1;
					if(r.numSites()>0 && r.numSites()>0){
						int z=r.topSite().score;
						for(int i=1; i<r.sites.size(); i++){
							SiteScore ss=r.sites.get(i);
							if(ss!=null && ss.score==z){x++;}
						}
					}
					if(r.ambiguous()){x=Tools.max(x, 2);}
				}
				if(x>=0){optionalTags.add("XM:i:"+x);}
			}

			//XS tag
			if(MAKE_XS_TAG){
				String xs=makeXSTag(r);
				if(xs!=null){
					optionalTags.add(xs);
					assert(r2==null || r.pairnum()!=r2.pairnum());
					//					assert(r2==null || !r2.mapped() || r.strand()==r2.strand() || makeXSTag(r2)==xs) :
					//						"XS problem:\n"+r+"\n"+r2+"\n"+xs+"\n"+makeXSTag(r2)+"\n";
				}
			}

			if(MAKE_MD_TAG){
				String md=makeMdTag(r.chrom, r.start, r.match, r.bases, scafloc, scaflen);
				if(md!=null){optionalTags.add(md);}
			}

			if(r.mapped() && MAKE_NH_TAG){
				if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && r.numSites()>1){
					optionalTags.add("NH:i:"+r.sites.size());
				}else{
					optionalTags.add("NH:i:1");
				}
			}

			if(MAKE_STOP_TAG && (perfect || (r.match!=null && r.bases!=null))){optionalTags.add(makeStopTag(pos, r.length(), cigar, perfect));}
			
			if(MAKE_LENGTH_TAG && (perfect || (r.match!=null && r.bases!=null))){optionalTags.add(makeLengthTag(pos, r.length(), cigar, perfect));}

			if(MAKE_IDENTITY_TAG && (perfect || r.match!=null)){optionalTags.add(makeIdentityTag(r.match, perfect));}

			if(MAKE_MATEQ_TAG && r.mateMapped()){
				optionalTags.add("YQ:i:"+toMapq(r.mate, null));
				optionalTags.add(makeIdentityTag(r.mate.match, r.mate.perfect()).replace('I', 'J'));
			}
			
			if(MAKE_SCORE_TAG && r.mapped()){optionalTags.add(makeScoreTag(r.mapScore));}

			if(MAKE_INSERT_TAG && r2!=null){
				if((r.mapped() && r.paired()) || r.originalSite!=null){
					optionalTags.add("X8:Z:"+r.insertSizeMapped(false)+(r.originalSite==null ? "" : ","+r.insertSizeOriginalSite()));
				}
			}
			if(MAKE_CORRECTNESS_TAG){
				final SiteScore ss0=r.originalSite;
				if(ss0!=null){
					optionalTags.add("X9:Z:"+(ss0.isCorrect(r.chrom, r.strand(), r.start, r.stop, 0) ? "T" : "F"));
				}
			}
		}

		if(READGROUP_ID!=null){
			assert(READGROUP_TAG!=null);
			optionalTags.add(READGROUP_TAG);
		}

		if(MAKE_CUSTOM_TAGS){
			int sites=r.numSites() + (r.originalSite==null ? 0 : 1);
			if(sites>0){
				ByteBuilder sb=new ByteBuilder();
				sb.append("X1:Z:");
				if(r.sites!=null){
					for(SiteScore ss : r.sites){
						sb.append('$');
						sb.append(ss.toText());
					}
				}
				if(r.originalSite!=null){
					sb.append('$');
					sb.append('*');
					sb.append(r.originalSite.toText());
				}
				optionalTags.add(sb.toString());
			}

			if(mapped){
				if(r.match!=null){
					byte[] match=r.match;
					if(!r.shortmatch()){
						match=Read.toShortMatchString(match);
					}
					optionalTags.add("X2:Z:"+new String(match));
				}

				optionalTags.add("X3:i:"+r.mapScore);
			}
			optionalTags.add("X5:Z:"+r.numericID);
			optionalTags.add("X6:i:"+(r.flags|(r.match==null ? 0 : Read.SHORTMATCHMASK)));
			if(r.copies>1){optionalTags.add("X7:i:"+r.copies);}
		}
		
		if(MAKE_TIME_TAG){
			assert(r.obj!=null && r.obj.getClass()==Long.class) : r.obj;
			optionalTags.add("X0:i:"+(r.obj==null ? 0 : r.obj));
		}
		
		if(MAKE_BOUNDS_TAG){
			String a=(r.mapped() ? inbounds ? "I" : "O" : "U");
			if(r2==null){
				optionalTags.add("XB:Z:"+a);
			}else{
				String b=(r2.mapped() ? inbounds2 ? "I" : "O" : "U");
				optionalTags.add("XB:Z:"+a+b);
			}
		}
		
		return optionalTags;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            ?            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Length of read bases */
	public int length(){
		assert((seq!=null && (seq.length!=1 || seq[0]!='*')) || cigar!=null) :
			"This program requires bases or a cigar string for every sam line.  Problem line:\n"+this+"\n";
		return seq==null ? calcCigarBases(cigar, true, false) : seq.length;
	}
	
//	public int length(boolean includeSoftClip){
//		assert((seq!=null && (seq.length!=1 || seq[0]!='*')) || cigar!=null) :
//			"This program requires bases or a cigar string for every sam line.  Problem line:\n"+this+"\n";
//		return seq==null ? calcCigarBases(cigar, includeSoftClip, false) : seq.length;
//	}
	
	public static int toMapq(Read r, SiteScore ss){
		assert(r!=null);
		int score=(ss==null ? r.mapScore : ss.slowScore);
		return toMapq(score, r.length(), r.mapped(), r.ambiguous());
	}
	
	public static int toMapq(int score, int length, boolean mapped, boolean ambig){
		if(!mapped || length<1){return 0;}
		
		if(ambig && PENALIZE_AMBIG){
			float max=3;
			float adjusted=(score*max)/(100f*length);
			return Tools.max(1, (int)Math.round(adjusted));
		}else{
			float score2=(score-length*40)*1.6f;
			float max=1.5f*((float)Tools.log2(length))+36;
			float adjusted=(score2*max)/(100f*length);
			return Tools.max(4, (int)Math.round(adjusted));
		}
	}
	
	
	public Read parseName(){
		try {
			String[] answer=qname.split("_");
			long id=Long.parseLong(answer[0]);
			int trueChrom=Gene.toChromosome(answer[1]);
			byte trueStrand=Byte.parseByte(answer[2]);
			int trueLoc=Integer.parseInt(answer[3]);
			int trueStop=Integer.parseInt(answer[4]);
//			for(int i=0; i<quals.length; i++){quals[i]-=33;}
//			Read r=new Read(seq.getBytes(), trueChrom, trueStrand, trueLoc, trueStop, qname, quals, false, id);
			Read r=new Read(seq, qual, qname, id, trueStrand, trueChrom, trueLoc, trueStop);
			return r;
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	
	public long parseNumericId(){
//		return Long.parseLong(qname.substring(0, qname.indexOf('_')));
		return Long.parseLong(qname.split("_")[1]);
	}
	
	public Read toRead(boolean parseCustom){
		return toRead(parseCustom, false);
	}
	
	public Read toRead(boolean parseCustom, boolean includeHardClip){
		
		SiteScore originalSite=null;
		long numericId_=0;
		boolean synthetic=false;
		
		if(parseCustom){
			
			
			CustomHeader h=new CustomHeader(qname, pairnum());
			
			numericId_=h.id;
			int trueChrom=h.bbchrom;
			byte trueStrand=(byte)h.strand;
			int trueLoc=h.bbstart;
			int trueStop=h.bbstop();
			
			originalSite=new SiteScore(trueChrom, trueStrand, trueLoc, trueStop, 0, 0);
			synthetic=true;
			
//			try {
//				String[] answer=qname.split("_");
//				numericId_=Long.parseLong(answer[0]);
//				int trueChrom=Gene.toChromosome(answer[1]);
//				byte trueStrand=Byte.parseByte(answer[2]);
//				int trueLoc=Integer.parseInt(answer[3]);
//				int trueStop=Integer.parseInt(answer[4]);
//				
//				originalSite=new SiteScore(trueChrom, trueStrand, trueLoc, trueStop, 0, 0);
//				synthetic=true;
//				
//			} catch (NumberFormatException e) {
//				System.err.println("Failed to parse "+qname);
//			} catch (NullPointerException e) {
//				System.err.println("Bad read with no name.");
//				return null;
//			}
		}
//		assert(false) : originalSite;
		
		
		if(Data.GENOME_BUILD>=0){
			
		}
		
		int chrom_=-1;
		byte strand_=strand();
		int start_=start(true, includeHardClip);
		int stop_=stop(start_, true, includeHardClip);
		assert(start_<=stop_) : start_+", "+stop_+"\n"+this+"\n";
		
		if(Data.GENOME_BUILD>=0){
			ScafLoc sc=null;
			if(RNAME_AS_BYTES){
				if(rname!=null && (rname.length!=1 || rname[0]!='*')){
					sc=Data.getScafLoc(rname);
					assert(sc!=null) : "Can't find scaffold in reference with name "+new String(rname)+"\n"+this;
				}
			}else{
				if(rnameS!=null && (rnameS.length()!=1 || rnameS.charAt(0)!='*')){
					sc=Data.getScafLoc(rnameS);
					assert(sc!=null) : "Can't find scaffold in reference with name "+new String(rnameS)+"\n"+this;
				}
			}
			if(sc!=null){
				chrom_=sc.chrom;
				start_+=sc.loc;
				stop_+=sc.loc;
			}
		}

////		byte[] quals=(qual==null || (qual.length()==1 && qual.charAt(0)=='*')) ? null : qual.getBytes();
////		byte[] quals=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual.clone();
//		byte[] quals=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual;
//		byte[] bases=seq==null ? null : seq.clone();
//		if(strand_==Gene.MINUS){//Minus-mapped SAM lines have bases and quals reversed
//			AminoAcid.reverseComplementBasesInPlace(bases);
//			Tools.reverseInPlace(quals);
//		}
//		Read r=new Read(bases, chrom_, strand_, start_, stop_, qname, quals, cs_, numericId_);
		
		final Read r;
		{
			byte[] seqX=(seq==null || (seq.length==1 && seq[0]=='*')) ? null : seq;
			byte[] qualX=(qual==null || (qual.length==1 && qual[0]=='*')) ? null : qual;
			String qnameX=(qname==null || qname.equals(stringstar)) ? null : qname;
			r=new Read(seqX, qualX, qnameX, numericId_, strand_, chrom_, start_, stop_);
		}
		
		r.setMapped(mapped());
		r.setSynthetic(synthetic);
//		r.setPairnum(pairnum()); //TODO:  Enable after fixing assertions that this will break in read input streams.
		if(originalSite!=null){
			r.originalSite=originalSite;
		}
		
		r.mapScore=mapq;
		r.setSecondary(!primary());
		
//		if(mapped()){
//			r.list=new ArrayList<SiteScore>(1);
//			r.list.add(new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0));
//		}
		
//		System.out.println(optional);
		if(optional!=null){
			for(String s : optional){
				if(s.equals("XT:A:R")){
					r.setAmbiguous(true);
				}else if(s.startsWith("X1:Z:")){
//					System.err.println("Found X1 tag!\t"+s);
					String[] split=s.split("\\$");
//					assert(false) : Arrays.toString(split);
					ArrayList<SiteScore> list=new ArrayList<SiteScore>(3);
					
					for(int i=1; i<split.length; i++){
//						System.err.println("Processing ss\t"+split[i]);
						String s2=split[i];
						SiteScore ss=SiteScore.fromText(s2);
						if(s2.charAt(0)=='*'){
							r.originalSite=ss;
						}else{
							list.add(ss);
						}
					}
//					System.err.println("List size = "+list.size());
					if(list.size()>0){r.sites=list;}
				}else if(s.startsWith("X2:Z:")){
					String s2=s.substring(5);
					r.match=s2.getBytes();
				}else if(s.startsWith("X3:i:")){
					String s2=s.substring(5);
//					r.mapScore=Integer.parseInt(s2); //Messes up generation of ROC curve
				}else if(s.startsWith("X5:Z:")){
					String s2=s.substring(5);
					r.numericID=Long.parseLong(s2);
				}else if(s.startsWith("X6:i:")){
					String s2=s.substring(5);
					r.flags=Integer.parseInt(s2);
				}else if(s.startsWith("X7:i:")){
					String s2=s.substring(5);
					r.copies=Integer.parseInt(s2);
				}else{
//					System.err.println("Unknown SAM field:"+s);
				}
			}
		}
//		assert(false) : CONVERT_CIGAR_TO_MATCH;
		if(r.match==null && cigar!=null && (CONVERT_CIGAR_TO_MATCH || cigar.indexOf('=')>=0)){
//			r.match=cigarToShortMatch(cigar, true);
			r.match=toShortMatch(false);
			
			if(r.match!=null){
				r.setShortMatch(true);
				if(Tools.indexOf(r.match, (byte)'B')>=0){
					boolean success=r.fixMatchB();
//					if(!success){r.match=null;}
//					assert(false) : new String(r.match);
				}
//				assert(false) : new String(r.match);
			}
//			assert(false) : new String(r.match);
//			System.err.println(">\n"+cigar+"\n"+(r.match==null ? "null" : new String(r.match)));
		}
//		assert(false) : new String(r.match);
		
//		System.err.println("Resulting read: "+r.toText());
		
		return r;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Aproximate length of result of SamLine.toText() */
	public int textLength(){
		int len=11; //11 tabs
		len+=(3+9+3+9);
		len+=(tlen>999 ? 9 : 3);
		
		len+=(qname==null ? 1 : qname.length());
		len+=rnameLen();
		len+=(rnext==null ? 1 : rnext.length);
		len+=(cigar==null ? 1 : cigar.length());
		len+=(seq==null ? 1 : seq.length);
		len+=(qual==null ? 1 : qual.length);
		
		if(optional!=null){
			len+=optional.size();
			for(String s : optional){len+=s.length();}
		}
		return len;
	}
	
	public ByteBuilder toText(){return toBytes((ByteBuilder)null);}
	
	public ByteBuilder toBytes(ByteBuilder bb){
		
		final int buflen=Tools.max(rnameLen(), (rnext==null ? 1 : rnext.length), (seq==null ? 1 : seq.length), (qual==null ? 1 : qual.length));
		
		if(bb==null){bb=new ByteBuilder(textLength()+4);}
		if(qname==null){bb.append('*').tab();}else{bb.append(qname).tab();}
		bb.append(flag).tab();
		if(RNAME_AS_BYTES){
			assert(!(rname==null && rnameS!=null));
			appendTo(bb, rname).tab();
		}else{
			assert(!(rname!=null && rnameS==null)) : RNAME_AS_BYTES+", "+rname+", "+rnameS;
			appendTo(bb, rnameS).tab();
		}
		bb.append(pos).tab();
		bb.append(mapq).tab();
		if(cigar==null){bb.append('*');}
		else{bb.append(cigar);}
		bb.tab();
		appendTo(bb, rnext).tab();
		bb.append(pnext).tab();
		bb.append(tlen).tab();
		
		if(mapped() && strand()==Shared.MINUS){
			appendReverseComplemented(bb, seq).tab();
			appendQualReversed(bb, qual);
		}else{
			appendTo(bb, seq).tab();
			appendQual(bb, qual);
		}
		
		if(optional!=null){
			for(String s : optional){
				bb.tab().append(s);
			}
		}
		return bb;
	}
	
	@Override
	public String toString(){return toBytes(null).toString();}
	
	private static ByteBuilder appendTo(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		return sb.append(a);
	}
	
	private static ByteBuilder appendTo(ByteBuilder sb, String a){
		if(a==null || a==stringstar || (a.length()==1 && a.charAt(0)=='*')){return sb.append('*');}
		return sb.append(a);
	}
	
	private static ByteBuilder appendReverseComplemented(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}

		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=AminoAcid.baseToComplementExtended[a[j]];}
		sb.length+=a.length;

		return sb;
	}
	
	private static ByteBuilder appendQual(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}

		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=0; j<a.length; i++, j++){buffer[i]=(byte)(a[j]+33);}
		sb.length+=a.length;

		return sb;
	}
	
	private static ByteBuilder appendQualReversed(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		
		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=(byte)(a[j]+33);}
		sb.length+=a.length;
		
		return sb;
	}
	
	/** Assumes a custom name including original location */
	public byte[] originalContig(){
//		assert(PARSE_CUSTOM);
		int loc=-1;
		int count=0;
		for(int i=0; i<qname.length() && loc==-1; i++){
			if(qname.charAt(i)=='_'){
				count++;
				if(count==6){loc=i;}
			}
		}
		if(loc==-1){
			return null;
		}
		return qname.substring(loc+1).getBytes();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------             Flag             ----------------*/
	/*--------------------------------------------------------------*/
	
//	Bit Description
//	0x1 template having multiple fragments in sequencing
//	0x2 each fragment properly aligned according to the aligner
//	0x4 fragment unmapped
//	0x8 next fragment in the template unmapped
//	0x10 SEQ being reverse complemented
//	0x20 SEQ of the next fragment in the template being reversed
//	0x40 the first fragment in the template
//	0x80 the last fragment in the template
//	0x100 secondary alignment
//	0x200 not passing quality controls
//	0x400 PCR or optical duplicate
//	0x800 supplementary alignment
	

	public static int makeFlag(Read r, Read r2, int fragNum, boolean sameScaf){
		int flag=0;
		if(r2!=null){
			flag|=0x1;

			if(r.mapped() && r.valid() && r.match!=null &&
					(r2==null || (sameScaf && r.paired() && r2.mapped() && r2.valid() && r2.match!=null))){flag|=0x2;}
			if(fragNum==0){flag|=0x40;}
			if(fragNum>0){flag|=0x80;}
		}
		if(!r.mapped()){flag|=0x4;}
		if(r2!=null && !r2.mapped()){flag|=0x8;}
		if(r.strand()==Shared.MINUS){flag|=0x10;}
		if(r2!=null && r2.strand()==Shared.MINUS){flag|=0x20;}
		if(r.secondary()){flag|=0x100;}
		if(r.discarded()){flag|=0x200;}
		return flag;
	}
	
	public boolean hasCigar() {
		return cigar!=null && cigar.length()>0 && cigar.charAt(0)!='*';
	}
	
	/** contains a cigar with X or = symbols */
	public boolean hasCigarXE() {
		if(cigar==null || cigar.length()<1 || cigar.charAt(0)=='*') {return false;}
		for(int i=0; i<cigar.length(); i++) {
			char c=cigar.charAt(i);
			if(c=='=' || c=='X') {return true;}
		}
		return false;
	}
	
	public boolean hasMate(){
		return (flag&0x1)==0x1;
	}
	
	public boolean properPair(){
		return (flag&0x2)==0x2;
	}
	
	public static boolean mapped(int flag){
		return (flag&0x4)!=0x4;
	}

	public static byte strand(int flag){
		return ((flag&0x10)==0x10 ? (byte)1 : (byte)0);
	}
	
	public boolean mapped(){
		return (flag&0x4)!=0x4;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}
	
	public boolean nextMapped(){
		return (flag&0x8)!=0x8;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}

	public byte strand(){
		return ((flag&0x10)==0x10 ? (byte)1 : (byte)0);
	}

	public byte mateStrand(){return nextStrand();}
	public byte nextStrand(){
		return ((flag&0x20)==0x20 ? (byte)1 : (byte)0);
	}
	
	public boolean firstFragment(){
		return (flag&0x40)==0x40;
	}
	
	public boolean lastFragment(){
		return (flag&0x80)==0x80;
	}
	
	public int pairnum(){
		return firstFragment() ? 0 : lastFragment() ? 1 : 0;
	}

	public boolean primary(){return (flag&0x100)==0;}
	public void setPrimary(boolean b){
		if(b){
			flag=flag&~0x100;
		}else{
			flag=flag|0x100;
		}
	}
	public void setMapped(boolean b){
		if(b){
			flag=flag&~0x4;
		}else{
			flag=flag|0x4;
		}
	}
	public void setFirstFragment(boolean b){
		if(b){
			flag=flag|0x40;
		}else{
			flag=flag&~0x40;
		}
	}
	public void setStrand(int strand){
		if(strand==1){
			flag=flag|0x10;
		}else{
			assert(strand==0);
			flag=flag&~0x10;
		}
	}
	
	public boolean discarded(){
		return (flag&0x200)==0x200;
	}
	
	public boolean duplicate(){
		return (flag&0x400)==0x400;
	}
	
	public boolean supplementary(){
		return (flag&0x800)==0x800;
	}
	
	public boolean leftmost(){
		if(!pairedOnSameChrom() || tlen==0){return true;}
		return tlen>0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             ?             ----------------*/
	/*--------------------------------------------------------------*/
	
//	/** Assumes rname is an integer. */
//	public int chrom(){
//		if(Data.GENOME_BUILD<0){return -1;}
//		HashMap sc
//	}
	
	public boolean ambiguous() {return mapped() && mapq<4;}
	
	/** Assumes rname is an integer. */
	public int chrom_old(){
		assert(false);
		if(!Tools.isDigit(rname[0]) && !Tools.isDigit(rname[rname.length-1])){
			if(warning){
				warning=false;
				System.err.println("Warning - sam lines need a chrom field.");
			}
			return -1;
		}
		assert(Shared.anomaly || '*'==rname[0] || (Tools.isDigit(rname[0]) && Tools.isDigit(rname[rname.length-1]))) :
			"This is no longer correct, considering that sam lines are named by scaffold.  They need a chrom field.\n"+new String(rname);
		if(rname==null || Arrays.equals(rname, bytestar) || !(Tools.isDigit(rname[0]) && Tools.isDigit(rname[rname.length-1]))){return -1;}
		//return Gene.toChromosome(new String(rname));
		//return Integer.parseInt(new String(rname)));
		final byte z='0';
		int x=rname[0]-z;
		for(int i=1; i<rname.length; i++){
			x=(x*10)+(rname[i]-z);
		}
		return x;
	}
	
	/** Returns the zero-based starting location of this read on the sequence. */
	public int start(boolean includeSoftClip, boolean includeHardClip){
		int x=countLeadingClip(cigar, includeSoftClip, includeHardClip);
		return pos-1-x;
	}
	
	/** Returns the zero-based stop location of this read on the sequence. */
	public int stop(int start, boolean includeSoftClip, boolean includeHardClip){
		if(!mapped() || cigar==null || cigar.charAt(0)=='*'){
//			return -1;
			return start+(seq==null ? 0 : Tools.max(0, seq.length-1));
		}
		int r=start+calcCigarLength(cigar, includeSoftClip, includeHardClip)-1;

//		assert(false) : start+", "+r+", "+calcCigarLength(cigar, includeHardClip);
//		System.err.println("start= "+start+", stop="+r);
		return r;
	}
	
	public int stop2(final int start, final boolean includeSoftClip, final boolean includeHardClip){
		if(mapped() && cigar!=null && cigar.charAt(0)!='*'){return stop(start, includeSoftClip, includeHardClip);}
//		return (seq==null ? -1 : start()+seq.length());
		return (seq==null ? -1 : start+seq.length);
	}
	
	public long numericId(){
		return 0;
	}
	
	/** This includes half-mapped pairs. */
	public boolean pairedOnSameChrom(){
//		assert(false) : (rname==null ? "nullX" : new String(rname))+", "+
//		(rnext==null ? "nullX" : new String(rnext))+", "+Tools.equals(rnext, byteequals)+", "+Arrays.equals(rname, rnext)+"\n"+this;
		if(RNAME_AS_BYTES){
			return Tools.equals(rnext, byteequals) || Arrays.equals(rname, rnext) || (/*pairnum()==1 &&*/ Tools.equals(rname, byteequals));
		}else{
			return Tools.equals(rnext, byteequals) || Tools.equals(rnameS, rnext) || (/*pairnum()==1 &&*/ stringequals.equals(rnameS));
		}
	}
	
	/** Assumes a custom name including original location */
	public int originalContigStart(){
//		assert(PARSE_CUSTOM);
		int loc=-1;
		int count=0;
		for(int i=0; i<qname.length() && loc==-1; i++){
			if(qname.charAt(i)=='_'){
				count++;
				if(count==5){loc=i;}
			}
		}
		if(loc==-1){
			return -1;
		}
		
		int sum=0;
		int mult=1;
		for(int i=loc+1; i<qname.length(); i++){
			char c=qname.charAt(i);
			if(!Tools.isDigit(c)){
				if(i==loc+1 && c=='-'){mult=-1;}
				else{break;}
			}else{
				sum=(sum*10)+(c-'0');
			}
		}
		return sum*mult;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	public int rnameLen(){
		return (rname==null ? rnameS==null ? 1 : rnameS.length() : rname.length);
	}

	public byte[] rname(){
		assert(RNAME_AS_BYTES);
		return rname;
	}
	public byte[] rnext(){return rnext;}

	public void setRname(byte[] x){assert(RNAME_AS_BYTES);rname=x;}
	public void setRnext(byte[] x){rnext=x;}

	public void setRname(String x){assert(!RNAME_AS_BYTES);rnameS=x;}
	public void setRnext(String x){rnext=(x==null ? null : x.getBytes());}
	
	public String rnameS(){return rnameS!=null ? rnameS : rname==null ? null : new String(rname);}
	public String rnextS(){return rnext==null ? null : new String(rnext);}
	
	public String rnamePrefix() {
		return (rnameS!=null ? toPrefix(rnameS) : toPrefix(rname));
	}
	
	private static String toPrefix(String s) {
		for(int i=0; i<s.length(); i++) {
			if(Character.isWhitespace(s.charAt(i))) {
				return s.substring(0, i);
			}
		}
		return s;
	}
	
	private static String toPrefix(byte[] s) {
		for(int i=0; i<s.length; i++) {
			if(Character.isWhitespace(s[i])) {
				return new String(s, 0, i);
			}
		}
		return new String(s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addOptionalTag(String s) {
		if(optional==null) {optional=new ArrayList<String>();}
		optional.add(s);
	}
	
	public String findTag(String prefix) {
		if(optional==null){return null;}
		for(String s : optional){
			if(s.startsWith(prefix)){return s;}
		}
		return null;
	}

	public String mdTag(){return findTag("MD:Z:");}
	public String mateqTag(){return findTag("YQ:i:");}
	public int parseIntFlag(String prefix) {
		String tag=findTag(prefix);
		return tag==null ? Integer.MIN_VALUE : Parse.parseInt(tag, 5);
	}
	public float parseFloatFlag(String prefix) {
		String tag=findTag(prefix);
		return tag==null ? -1 : Parse.parseFloat(tag, 5);
	}
	public int mateq() {return parseIntFlag("YQ:i:");}
	public float mateID() {return parseFloatFlag("YJ:f:");}

	public int setScafnum(ScafMap scafMap) {
		assert(scafnum<0);
		
		String name=null;
		if(mapped() || (rname!=null && rname!=byteequals && rname!=bytestar)){
			name=rnameS();
		}else if(nextMapped() && rnext!=null && rnext!=byteequals && rnext!=bytestar){
			name=new String(rnext);
		}
		if(name!=null){scafnum=scafMap.getNumber(name);}
		return scafnum;
	}
	
	public long countBytes(){
		long sum=76;
		sum+=(cigar==null ? 0 : cigar.length()*2+16);
		sum+=(optional==null ? 0 : optional.size()*32+16);
		sum+=(rname==null ? 0 : rname.length+16);
		sum+=(rnext==null ? 0 : rnext.length+16);
		return sum;
	}
	
	public String qname;
	public int flag;
	public int pos;
	public int mapq;
	public String cigar;
	public int pnext;
	public int tlen;
	public byte[] seq;
	public byte[] qual;
	public ArrayList<String> optional;
	public byte[] mdTag;
	
	public Object obj;
	public int scafnum=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private byte[] rname;
	private byte[] rnext;
	
	private String rnameS;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final String stringstar="*";
	private static final String stringequals="=";
	private static final byte[] bytestar=new byte[] {(byte)'*'};
	private static final byte[] byteequals=new byte[] {(byte)'='};
	private static final String XSPLUS="XS:A:+", XSMINUS="XS:A:-";
//	private static final double inv100=0.01d;
//	private static float minratio=0.4f;

	private static boolean warning=System.getProperty("user.dir").contains("/bushnell/");
	
	/*--------------------------------------------------------------*/
	/*----------------     Public Static Fields     ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean makeReadgroupTags(){
		return READGROUP_ID!=null || READGROUP_CN!=null || READGROUP_DS!=null || READGROUP_DT!=null ||
				READGROUP_FO!=null || READGROUP_KS!=null || READGROUP_LB!=null || READGROUP_PG!=null || 
				READGROUP_PI!=null || READGROUP_PL!=null || READGROUP_PU!=null || READGROUP_SM!=null ||
				READGROUP_TAG!=null;
	}
	
	public static boolean makeOtherTags(){
		if(NO_TAGS){return false;}
		return MAKE_AM_TAG || MAKE_NM_TAG || MAKE_SM_TAG || MAKE_XM_TAG || MAKE_XS_TAG || MAKE_AS_TAG ||
				MAKE_NH_TAG || MAKE_TOPHAT_TAGS || MAKE_IDENTITY_TAG || MAKE_SCORE_TAG || MAKE_STOP_TAG || MAKE_LENGTH_TAG ||
				MAKE_CUSTOM_TAGS || MAKE_INSERT_TAG || MAKE_CORRECTNESS_TAG || MAKE_TIME_TAG || MAKE_BOUNDS_TAG || MAKE_MATEQ_TAG;
	}
	
	public static boolean makeAnyTags(){
		return makeReadgroupTags() || makeOtherTags();
	}
	
	public static String READGROUP_ID=null;
	public static String READGROUP_CN=null;
	public static String READGROUP_DS=null;
	public static String READGROUP_DT=null;
	public static String READGROUP_FO=null;
	public static String READGROUP_KS=null;
	public static String READGROUP_LB=null;
	public static String READGROUP_PG=null;
	public static String READGROUP_PI=null;
	public static String READGROUP_PL=null;
	public static String READGROUP_PU=null;
	public static String READGROUP_SM=null;
	
	public static String READGROUP_TAG=null;
	
	/** Turn this off for RNAseq or long indels */
	public static boolean MAKE_MD_TAG=false;
	
	public static boolean NO_TAGS=false;
	
	public static boolean MAKE_AM_TAG=true;
	public static boolean MAKE_NM_TAG=true;
	public static boolean MAKE_SM_TAG=false;
	public static boolean MAKE_XM_TAG=false;
	public static boolean MAKE_XS_TAG=false;
	public static boolean MAKE_XT_TAG=true;
	public static boolean MAKE_AS_TAG=false; //TODO: Alignment score from aligner
	public static boolean MAKE_NH_TAG=false;
	public static boolean MAKE_TOPHAT_TAGS=false;
	public static boolean XS_SECONDSTRAND=false;
	public static boolean MAKE_IDENTITY_TAG=false;
	public static boolean MAKE_SCORE_TAG=false;
	public static boolean MAKE_STOP_TAG=false;
	public static boolean MAKE_LENGTH_TAG=false;
	public static boolean MAKE_CUSTOM_TAGS=false;
	public static boolean MAKE_INSERT_TAG=false;
	public static boolean MAKE_CORRECTNESS_TAG=false;
	public static boolean MAKE_TIME_TAG=false;
	public static boolean MAKE_BOUNDS_TAG=false;
	public static boolean MAKE_MATEQ_TAG=false;
	
	public static boolean PENALIZE_AMBIG=true;
	public static boolean CONVERT_CIGAR_TO_MATCH=true;
	public static boolean SOFT_CLIP=true;
	public static boolean SECONDARY_ALIGNMENT_ASTERISKS=true;
	/** OK to use the "setFrom" function which uses the old SamLine instead of translating the read, if a genome is not loaded. */
	public static boolean SET_FROM_OK=false;
	/** For paired reads, keep original names rather than changing read2's name to match read1 */
	public static boolean KEEP_NAMES=false;
	public static float VERSION=1.4f;
	/** Tells program when to use 'N' rather than 'D' in cigar strings */
	public static int INTRON_LIMIT=Integer.MAX_VALUE;
	public static boolean RNAME_AS_BYTES=true;//Effect on speed is negligible for pileup...
	
	/** Prefer MD tag over reference for translating cigar strings to match */
	public static boolean PREFER_MDTAG=false;
	/** Determine whether cigar X means match N or S.
	 * This makes sam loading substantially slower. */
	public static boolean FIX_MATCH_NS=false;
	
	public static boolean setxs=false;
	public static boolean setintron=false;
	
//	/** SSAHA2 incorrectly calculates the start position of reads with soft-clipped starts, and needs this enabled. */
//	public static boolean SUBTRACT_LEADING_SOFT_CLIP=true;
	/** Sort header scaffolds in alphabetical order to be more compatible with Tophat */
	public static boolean SORT_SCAFFOLDS=false;

	/** qname */
	public static boolean PARSE_0=true;
	/** rname */
	public static boolean PARSE_2=true;
	/** cigar */
	public static boolean PARSE_5=true;
	/** rnext */
	public static boolean PARSE_6=true;
	/** pnext */
	public static boolean PARSE_7=true;
	/** tlen */
	public static boolean PARSE_8=true;
	/** qual */
	public static boolean PARSE_10=true;
	public static boolean PARSE_OPTIONAL=true;
	public static boolean PARSE_OPTIONAL_MD_ONLY=false;
	public static boolean PARSE_OPTIONAL_MATEQ_ONLY=false;
	
	public static boolean FLIP_ON_LOAD=true;
	
	public static boolean verbose=false;
	
}
