package stream;

import java.io.Serializable;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.Gene;
import dna.ScafLoc;
import parse.LineParser1;
import parse.Parse;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import simd.Vector;
import structures.ByteBuilder;
import var2.ScafMap;
import var2.Scaffold;


/**
 * Mutable SAM alignment fields with Read conversion, CIGAR helpers and text output.
 * Parsing and output follow process-wide flags; this is not a complete SAM validator.
 * qual stores numeric Phred scores, with ASCII-33 conversion at text boundaries.
 * Parsing may reverse-complement sequence and reverse qualities under FLIP_ON_LOAD;
 * setters retain the orientation supplied by the caller. Configure name-storage and
 * orientation policies consistently before constructing and serializing records.
 * Copies share arrays, optional tags and auxiliary objects. Callers coordinate mutation,
 * borrowed references and changes to global policy when records are shared.
 *
 * @author Brian Bushnell
 */
public class SamLine implements Serializable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates default-valued fields with scafnum=-1; does not create a validated alignment. */
	public SamLine(){}

	/** Shallow-copies record state, sharing arrays, optional tags and auxiliary objects.
	 * @param sl Nonnull source; no normalization or independent copies are made */
	public SamLine(SamLine sl){setFrom(sl);}

	/** Builds fields from a Read using current name, alignment and optional-tag policies.
	 * With no Data.scaffoldLocs and a retained r1.samline, asserts SET_FROM_OK and shallow-copies
	 * that record, returning before rebuilding fields or applying fragNum. Otherwise mapped
	 * coordinates depend on loaded Data scaffold metadata. Multi-scaffold alignments can clear
	 * mapped/paired flags and match arrays on r1 or its mate. Sequence and numeric qualities
	 * normally share r1's arrays; configured secondary-record suppression omits them.
	 * Trailing-clip and out-of-bounds indel helpers scan original match arrays and expect
	 * expanded operations in the scanned regions. Later short-match expansion is local to
	 * CIGAR generation; optional-tag generation also receives the original Read objects.
	 * @param r1 Nonnull read; the rebuilding path requires an identifier and compatible mapping metadata
	 * @param fragNum Zero-based mate selector for rebuilt flags, conventionally 0 or 1 */
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
				//TODO [STR-334]: this uses r1's scaffold length, so it can scan extra mate suffix
				//operations on different scaffolds or when r1 is unmapped. The result affects only
				//pos1_mate, used for sameScaf TLEN/assertions where lengths agree; no output discrepancy is established.
				int tclippedIndels=countTrailingIndels(b2, scaflen, r2.match);
				if(verbose){
					System.err.println("leadingClip="+clip);
					System.err.println("clippedDels="+clippedIndels);
				}
				pos0_mate=(a2+1)+clip+clippedIndels;
				pos1_mate=(b2+1)-tclip-tclippedIndels;
				if(pos1_mate>scaflen2){pos1_mate=scaflen2;}//[stream/SamLine#001] FIXED: was 'if(pos1_mate>scaflen){pos1=scaflen;}' - clamped the WRONG variable (pos1, r1's stop) on the wrong bound (scaflen vs scaflen2); copy-paste of the r1 clamp above. pos1_mate feeds TLEN only when sameScaf (where scaflen==scaflen2).

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

		optional=makeOptionalTags(r1, r2, perfect, scafloc, scaflen, inbounds, inbounds2);
//		assert(r.pairnum()==1) : "\n"+r.toText(false)+"\n"+this+"\n"+r2;
	}

	/** Parses an alignment under PARSE_* selection flags, mutating the parser's field bounds.
	 * FLAG, POS, MAPQ and SEQ are read unconditionally; selected arrays are copied from the
	 * parser buffer and strings decoded as US-ASCII. Text QUAL is converted to numeric Phred.
	 * Mapped reverse-strand records are flipped in place when FLIP_ON_LOAD is enabled.
	 * The optional MD-prefix filter takes precedence over the YQ-prefix filter. Name trimming
	 * then follows TRIM_READ_DESCRIPTION; this constructor does not completely validate SAM.
	 * @param lp Nonnull parser containing a tab-delimited alignment, not a header line */
	public SamLine(LineParser1 lp){
		assert(!lp.startsWith('@')) : "Tried to make a SamLine from a header: "+lp.toString();

		if(PARSE_0){qname=lp.parseString(0);}
		flag=lp.parseInt(1);
		if(PARSE_2){
			lp.setBounds(2);
			boolean isStar=lp.currentTermEquals(star);
			if(RNAME_AS_BYTES){
				rname=(isStar) ? null : lp.parseByteArray(2);
			}else{
				rnameS=(isStar) ? null : lp.parseString(2);
			}
		}
		pos=lp.parseInt(3);
		mapq=lp.parseInt(4);
		if(PARSE_5){
			int len=lp.setBounds(5);
			setCigar(lp.currentTermEquals(star) ? null : lp.parseStringFromCurrentField());
		}
		if(PARSE_6){
			int len=lp.setBounds(6);
			rnext=lp.currentTermEquals(star) ? null : lp.parseByteArrayFromCurrentField();
		}

		if(PARSE_7){
			int len=lp.setBounds(7);
			pnext=(len<1 || lp.currentTermEquals(star)) ? 0 : lp.parseIntFromCurrentField();
		}
		if(PARSE_8){tlen=lp.parseInt(8);}
		{
			int len=lp.setBounds(9);
			seq=(len<1 || lp.currentTermEquals(star)) ? null : lp.parseByteArrayFromCurrentField();
		}
		if(PARSE_10){
			int len=lp.setBounds(10);
			qual=(len<1 || lp.currentTermEquals(star)) ? null : lp.parseByteArrayFromCurrentField();
		}

		assert((seq==bytestar)==(Tools.equals(seq, bytestar)));
		assert((qual==bytestar)==(Tools.equals(qual, bytestar)));

		if(FLIP_ON_LOAD && mapped() && strand()==Shared.MINUS){
			if(seq!=null && qual!=bytestar){Vector.reverseComplementInPlaceFast(seq);}
			if(qual!=null && qual!=bytestar){Vector.reverseInPlace(qual);}
		}

		if(qual!=null && qual!=bytestar){
//			for(int i=0; i<qual.length; i++){qual[i]-=33;}
			Vector.add(qual, (byte)(-33));
		}

		if(PARSE_OPTIONAL && lp.terms()>11){
			if(PARSE_OPTIONAL_MD_ONLY){
				optional=new ArrayList<String>(1);
				for(int i=11, terms=lp.terms(); i<terms; i++){
					if(lp.termStartsWith("MD:", i)){
						String s=lp.parseString(i);
						optional.add(s);
//						lp.incrementA(5);
//						mdTag=lp.parseByteArrayFromCurrentField();//Not really needed
					}
				}
			}else if(PARSE_OPTIONAL_MATEQ_ONLY){
				optional=new ArrayList<String>(1);
				for(int i=11, terms=lp.terms(); i<terms; i++){
					if(lp.termStartsWith("YQ:", i)){
						String s=lp.parseString(i);
						optional.add(s);
					}
				}
			}else{
				optional=new ArrayList<String>(lp.terms()-11);
				for(int i=11, terms=lp.terms(); i<terms; i++){
					String s=lp.parseString(i);
					optional.add(s);
				}
			}
		}

		trimNames();
	}

	/*--------------------------------------------------------------*/
	/*-------------------        Methods        --------------------*/
	/*--------------------------------------------------------------*/

	/** Copies stored scalar fields and references, including both name representations and cached metadata.
	 * No arrays, optional-tag lists or payload objects are cloned. Requires nonnull sl. */
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
		//[stream/SamLine#004] FIXED: previously omitted mdTag/obj/scafnum (later-added-field-not-copied pattern). A copy must reproduce ALL state.
		mdTag=sl.mdTag;
		obj=sl.obj;
		scafnum=sl.scafnum;
	}

	/** Applies input-side TRIM_READ_DESCRIPTION to QNAME and the active RNAME/RNEXT fields.
	 * Tools allocates a prefix array when Character whitespace is found.
	 * Name setters also canonicalize missing/equal text tokens; output TRIM_* flags are separate. */
	public void trimNames(){
		if(Shared.TRIM_READ_DESCRIPTION){
			if(RNAME_AS_BYTES){
				setRname(Tools.trimToWhitespace(rname()));
				setRnext(Tools.trimToWhitespace(rnext()));
			}else{
				setRnameS(Tools.trimToWhitespace(rnameS()));
				setRnext(Tools.trimToWhitespace(rnext()));
			}
			qname=(Tools.trimToWhitespace(qname));
		}
	}

	/** Tests only the character alphabet: M, = and Character digits, with nonempty input.
	 * This does not validate alternating positive counts and operations; null/empty returns false. */
	public boolean cigarContainsOnlyME(){
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

	/**
	 * Returns this record's reference span plus explicitly requested clipping contributions.
	 * @param includeSoftClip Whether to include soft-clipped bases
	 * @param includeHardClip Whether to include hard-clipped bases
	 * @return Reference count from calcCigarLength, or zero for null CIGAR
	 */
	public int calcCigarLength(boolean includeSoftClip, boolean includeHardClip){return calcCigarLength(cigar, includeSoftClip, includeHardClip);}

	/**
	 * Returns this record's query count with explicitly requested clipping contributions.
	 * @param includeSoftClip Whether to include soft-clipped bases
	 * @param includeHardClip Whether to include hard-clipped bases
	 * @return Query count from calcCigarReadLength, or zero for null CIGAR
	 */
	public int calcCigarReadLength(boolean includeSoftClip, boolean includeHardClip){return calcCigarReadLength(cigar, includeSoftClip, includeHardClip);}

	/** Sums M/=/X/I query counts for a mapped record, excluding clips and D/N gaps.
	 * Returns zero for unmapped records or null CIGAR. P/unknown operations and a trailing
	 * digit run contribute nothing; this counter does not validate CIGAR grammar. */
	public int mappedNonClippedBases(){
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

	/** Tests whether the default ScafMap resolves this record's reference name to nonnull bases.
	 * Lookup retains ScafMap's missing/ambiguous-name diagnostics. Does not check query
	 * sequence availability, alignment bounds or mapping flags. */
	private boolean refLoadedForThisRead(){
		final ScafMap sm=ScafMap.defaultScafMap();
		if(sm==null){return false;}
		final Scaffold scaf=sm.getScaffold(rnameS());
		return scaf!=null && scaf.bases!=null;
	}

	/** Estimates identity under this class's CIGAR/MD/reference policy; requires nonnull CIGAR.
	 * The direct path uses =/(=+X+D+I), with denominator at least one; N and clips are omitted.
	 * M alongside =/X is omitted as a BBMap no-call marker. M without =/X first tries
	 * toShortMatch(false) when MD or reference bases are available, using m/(m+S+I+D)
	 * from the returned match array. That path has toShortMatch's sequence-mutation contract.
	 * If no match array is obtained, asserts M_CIGARS_OK through the terminating diagnostic;
	 * permissive execution counts M as matches. This is not a general CIGAR validator.
	 * @return Identity fraction for consistent, non-overflowing counts; zero for a zero denominator */
	public float calcIdentity(){
		assert(cigar!=null);
		int match=0, other=0, mCount=0;
		boolean foundM=false, foundEX=false;

		int current=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='='){
					match+=current; foundEX=true;
				}else if(c=='M'){
					mCount+=current; foundM=true; //Resolved below: N-base (mixed with =/X) or ambiguous (M-only)
				}else if(c=='X'){
					other+=current; foundEX=true;
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
		//M alongside =/X is BBMap's N-base marker (neither match nor mismatch) -> correctly counted as nothing above.
		//M with NO =/X is an ambiguous match/mismatch cigar (e.g. minimap2): the cigar alone can't split match from sub.
		if(foundM && !foundEX){
			//Ambiguous M-only cigar (e.g. minimap2). Resolve the match/sub split via the MD tag or the loaded
			//reference -- exactly what toShortMatch does. Crash loud only if NEITHER is available and the tool
			//did not opt into permissive M handling.
			if(mdTag()!=null || refLoadedForThisRead()){
				final byte[] sm=toShortMatch(false); //resolves M via MD tag or ScafMap reference
				if(sm!=null){
					final int[] c=Read.countMatchSymbols(sm); //{m,S,C,N,I,D}
					return c[0]/(float)Tools.max(c[0]+c[1]+c[4]+c[5], 1); //m/(m+S+I+D)
				}
			}
			assert(M_CIGARS_OK) : KillSwitch.assertDie("Cannot compute identity for an ambiguous M-only cigar (no '='/'X'), "
					+ "and no MD tag or loaded reference to resolve it. Provide =/X cigars or an MD tag, load the reference, "
					+ "or use a tool that permits M (sets SamLine.M_CIGARS_OK). cigar="+cigar);
			match+=mCount; //M_CIGARS_OK (e.g. QuickBin): coverage-focused, count aligned M as match so the read is not dropped.
		}
		return match/(float)Tools.max(match+other, 1);
	}

	/** Counts X operations, resolving M-only alignments through toShortMatch(false) when possible.
	 * M mixed with =/X adds no substitutions. MD/reference resolution may temporarily modify
	 * seq as documented by toShortMatch. With no resulting match array, an M-only CIGAR
	 * requires M_CIGARS_OK under assertions; permissive execution returns zero substitutions.
	 * @return Substitution count, or zero for null CIGAR; does not validate CIGAR grammar */
	public int countSubs(){
		if(cigar==null){return 0;}

		int current=0;
		int subs=0;
		boolean foundM=false, foundEX=false;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(c=='X'){
					subs+=current; foundEX=true;
				}else if(c=='='){
					foundEX=true;
				}else if(c=='M'){
					foundM=true;
				}
				current=0;
			}
		}
		//Ambiguous M-only cigar (M with no =/X): substitutions are uncomputable from the cigar (see calcIdentity).
		if(foundM && !foundEX){
			if(mdTag()!=null || refLoadedForThisRead()){//Fixable via MD tag or loaded reference.
				final byte[] sm=toShortMatch(false);
				if(sm!=null){return Read.countMatchSymbols(sm)[1];}//S = substitutions
			}
			assert(M_CIGARS_OK) : KillSwitch.assertDie("Cannot count substitutions for an ambiguous M-only cigar (no '='/'X'), "
					+ "and no MD tag or loaded reference to resolve it. Provide =/X cigars or an MD tag, load the reference, "
					+ "or use a tool that permits M (sets SamLine.M_CIGARS_OK). cigar="+cigar);
			//M_CIGARS_OK and unresolvable: permissive, report 0 subs (subs is already 0 here; coverage-focused).
		}
		return subs; //Mixed M+=/X: M contributes no subs (N-bases); X count stands.
	}

	/** Converts CIGAR to newly allocated short BBTools match text, optionally resolving ambiguous M.
	 * Maps =/X/D-or-N/I/S to m/S/D/I/C and omits hard clips. With allowM=false, M mixed
	 * with =/X becomes N; M alone is resolved using MD or the default ScafMap reference.
	 * An MD tag with zero substitutions supplies the fast path for M-only input. Otherwise
	 * PREFER_MDTAG selects MD when present; usable reference/query bases are preferred by
	 * default, with MD fallback. FIX_MATCH_NS may also reclassify X no-calls when no M exists.
	 * Correction may reverse-complement stored seq on entry and again on normal completion
	 * for a minus-strand record with reference/query bases, independently of FLIP_ON_LOAD.
	 * Callers must supply the expected sequence orientation and coordinate shared-array access.
	 * Missing M-resolution context uses existing terminating assertions; absent query bases
	 * without usable MD can leave unresolved match symbols. This is not complete SAM validation.
	 * @param allowM If true, delegates immediately to cigarToShortMatch_old, encoding M as N
	 * without MD/reference correction or sequence mutation
	 * @return Short match array; null for absent/'*' CIGAR and, on the resolving path, P or
	 * no counted operations, with additional null outcomes in unsupported resolution cases */
	public final byte[] toShortMatch(boolean allowM){
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
					assert(!fixMatchSubs || mdTag!=null || map!=null) : KillSwitch.assertDie("Encountered a read with 'M' in cigar string but no MD tag and no ScafMap loaded.\n"
							+ "This can normally be resolved by adding the flag ref=file, where file is the fasta file to which the reads were mapped.\n\n"+this);
					if(map==null){
						refBases=null;
					}else{
						Scaffold scaf=map.getScaffold(rnameS());
						assert(!fixMatchSubs || mdTag!=null || scaf!=null) : KillSwitch.assertDie("Encountered a read with 'M' in cigar string but no scaffold loaded for "+rnameS());
						if(scaf==null){
							refBases=null;
						}else{
							refBases=scaf.bases;
							assert(!fixMatchSubs || mdTag!=null || refBases!=null) : KillSwitch.assertDie("Encountered a read with 'M' in cigar string but no ref bases loaded for "+rnameS());
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
		if(refBases!=null && seq!=null){
//			bases=(strand()==1) ? AminoAcid.reverseComplementBases(seq) : seq;//Why not reverse in place?
			if(strand()==1){Vector.reverseComplementInPlaceFast(seq);}
			bases=seq;
		}else{bases=null;}

//		System.err.println("Block 5.");//123
		final byte[] longmatch=Read.toLongMatchString(match0);

//		final int noCalls=((foundE || foundX) && foundM ? 1 : seq==null ? -1 : Read.countNocalls(seq));

//		System.err.println("Block 6");//123

		if(mdTag!=null && (refBases==null || PREFER_MDTAG || bases==null)){
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
		}else if(refBases!=null && bases!=null){
			final int refStart=start(true, false);
			fixMatch(bases, refBases, longmatch, refStart, false);
		}

//		else if(refBases!=null){
//			final int refStart=start(true, false);
//			fixMatch(bases, refBases, longmatch, refStart, false);
//		}
		else{
			//Null bases (e.g. SEQ=* secondary/supplementary) with no usable MD tag: can't compute =/X; leave longmatch as-is.
		}

		final byte[] match=Read.toShortMatchString(longmatch);

		if(bases!=null && strand()==1){Vector.reverseComplementInPlace(seq);}

//		System.err.println("Block 7.");//123
//		System.err.println("Returning "+new String(match));//123
		return match;
	}

	/** Returns the first optional string starting with X2:Z:, including its prefix.
	 * Null list/no match returns null; does not decode or validate match text. */
	public String matchTag(){
		if(optional==null){return null;}
		for(String s : optional){
			if(s.startsWith("X2:Z:")){
				return s;
			}
		}
		return null;
	}

	/**
	 * Selects a shared XS strand tag only when r is mapped and this CIGAR contains N.
	 * Starts from r's strand, inverts for a nonzero pair number, then for XS_SECONDSTRAND.
	 * This encodes the configured library convention, not an inferred transcript strand.
	 * @param r Nonnull read supplying mapping status, strand and pair number
	 * @return Shared XSPLUS/XSMINUS string, or null when the mapping/CIGAR condition fails
	 */
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

	/**
	 * Builds a new optional-tag list from r, this record's CIGAR/POS/MAPQ and global policies.
	 * Does not assign this.optional, copy existing optional tags, or deduplicate generated tags.
	 * NO_TAGS returns null. Unmapped r also returns null unless read-group, custom or timing
	 * output is enabled; bounds output alone does not bypass that early return.
	 * Mapping-dependent tags are generated only for mapped r. Tophat tags take precedence over
	 * the XM branch. NM and MD generation expect expanded r.match and consistent query/reference
	 * metadata. Custom X2 output separately honors r.shortmatch(). Selected tags consult both
	 * r.mate and the supplied r2; callers must keep those references consistent.
	 * Read-group output requires READGROUP_TAG and timing output asserts a Long r.obj.
	 * @param r Nonnull source read, including secondary/unmapped reads when supported by selected tags
	 * @param r2 Mate context, possibly null
	 * @param perfect Caller-supplied perfect status, used by NM and custom tag shortcuts
	 * @param scafloc Zero-based chromosome scaffold start for MD generation
	 * @param scaflen Scaffold length for MD generation
	 * @param inbounds Caller-supplied read bounds status for XB
	 * @param inbounds2 Caller-supplied mate bounds status for XB
	 * @return New list (possibly empty), or null at the global/unmapped early returns
	 */
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

						if(b=='D'){delsCurrent++;}else{
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
				if(perfect){optionalTags.add("NM:i:0");}else if(r.match!=null){optionalTags.add("NM:i:"+(nm));}
			}
			if(MAKE_SM_TAG){optionalTags.add("SM:i:"+mapq);}
			if(MAKE_AM_TAG){optionalTags.add("AM:i:"+Data.min(mapq, r2==null ? mapq : (r2.mapped() ? Data.max(1, r2.mapScore/r2.length()) : 0)));}
			if(MAKE_DUAL_MAPQ_TAGS && r.primary()){
				final int loose=NeuralMapqCache.getLoose(r), strict=NeuralMapqCache.getStrict(r);
				assert((loose<0)==(strict<0)) : KillSwitch.assertDie("QL/QS must be present together; loose="+loose+", strict="+strict+", read="+r.id);
				if(loose>=0){optionalTags.add("QL:i:"+loose); optionalTags.add("QS:i:"+strict);}
			}

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
					optionalTags.add("X2:Z:"+new String(match, StandardCharsets.US_ASCII));
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


	/** Returns stored sequence length, otherwise query CIGAR count including S but excluding H.
	 * Asserts that sequence is present and not '*' or that CIGAR exists; does not validate
	 * sequence/CIGAR agreement. Empty stored arrays return zero. */
	public int length(){
		assert((seq!=null && (seq.length!=1 || seq[0]!='*')) || cigar!=null) :
			"This program requires bases or a cigar string for every sam line.  Problem line:\n"+this+"\n";
		return seq==null ? calcCigarBases(cigar, true, false) : seq.length;
	}

	/** Uses stored sequence length, otherwise query CIGAR count including S/excluding H, otherwise zero. */
	public int lengthOrZero(){return seq!=null ? seq.length : cigar!=null ? calcCigarBases(cigar, true, false) : 0;}

	/** Returns a rough capacity estimate from fixed overhead, sequence, QNAME and CIGAR.
	 * Requires nonnull QNAME; omits several fields and tags, so this is not exact BAM size. */
	public int estimateBamLength(){return 40+(seq==null ? 1 : seq.length)+qname.length()+(cigar==null ? 1 : cigar.length()*2);}


	/** Parses the legacy underscore-separated ID/chromosome/strand/start/stop name fields.
	 * Shares seq/qual with the new Read; configured Read construction may modify those arrays.
	 * Only NumberFormatException is caught, printed and converted to null. Missing/null fields
	 * and other construction failures are not covered by that fallback; no CustomHeader is used.
	 * @return New Read carrying the parsed coordinates, or null after a caught numeric parse failure */
	public Read parseName(){
		try{
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
		}catch(NumberFormatException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}

	/** Parses underscore field index one as a long; unlike parseName, this does not use field zero.
	 * Requires a nonnull name and parseable field; failures propagate. */
	public long parseNumericId(){
//		return Long.parseLong(qname.substring(0, qname.indexOf('_')));
		return Long.parseLong(qname.split("_")[1]);
	}

	/**
	 * Delegates to toRead(parseCustom, false), excluding hard clips from coordinates.
	 * @param parseCustom Whether to parse custom BBTools naming format
	 * @return Read object representation
	 */
	public Read toRead(boolean parseCustom){return toRead(parseCustom, false);}

	/**
	 * Converts this SamLine to a Read object with detailed options.
	 * Uses stored sequence orientation, coordinate helpers, and optional tag parsing.
	 * Retains sequence and numeric-quality arrays; configured Read validation may modify them.
	 * A loaded Data genome can translate reference-local coordinates into chromosome coordinates.
	 * Recognized custom tags are applied in list order, even when parseCustom is false; X6 replaces
	 * the Read flags. No mate link or SamLine back-reference is attached, and the normal SAM pair
	 * number assignment is disabled. MAPQ becomes mapScore; selected optional tags may replace
	 * numeric ID, flags, copies, sites or match. Missing match may be generated from CIGAR under
	 * CONVERT_CIGAR_TO_MATCH or when '=' occurs, with toShortMatch's mutation/dependency contract.
	 * @param parseCustom Whether to parse custom BBTools naming format
	 * @param includeHardClip Whether to include hard-clipped bases in coordinates
	 * @return Read object with alignment and metadata
	 */
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
//			Vector.reverseComplementInPlace(bases);
//			Tools.reverseInPlace(quals);
//		}
//		Read r=new Read(bases, chrom_, strand_, start_, stop_, qname, quals, cs_, numericId_);

		final Read r;
		{
			byte[] seqX=(seq==null || (seq.length==1 && seq[0]=='*')) ? null : seq;
			//Numeric Phred 42 is a score; only the legacy missing-array identity is a sentinel.
			byte[] qualX=(qual==null || qual==bytestar) ? null : qual;
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
		r.setSecondary(!nonSecondary());
		r.setSupplementary(supplementary());

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


	/** Returns a historical SAM text capacity estimate, not an exact serialized length.
	 * Uses fixed numeric widths and does not account precisely for output trimming or all values. */
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

	/** Returns a new builder containing SAM alignment text without a trailing newline. */
	public ByteBuilder toText(){return toBytes((ByteBuilder)null);}

	/** Appends SAM alignment text under current output-name and orientation policies.
	 * Appends to existing builder contents and omits a final newline. TRIM_QNAME/TRIM_RNAME
	 * affect emitted names without replacing stored fields. Mapped minus-strand sequence
	 * and qualities are emitted in reverse order only when FLIP_ON_LOAD is enabled; sequence
	 * is complemented and numeric qualities gain 33. Record arrays are not reversed in place.
	 * The active RNAME_AS_BYTES representation must be consistent with the stored name.
	 * @param bb Caller-owned destination, allocated when null; storage must not alias record arrays
	 * @return Destination builder with the record appended */
	public ByteBuilder toBytes(ByteBuilder bb){

		final int buflen=Tools.max(rnameLen(), (rnext==null ? 1 : rnext.length), (seq==null ? 1 : seq.length), (qual==null ? 1 : qual.length));

		if(bb==null){bb=new ByteBuilder(textLength()+4);}
		if(qname==null){bb.append('*').tab();}else{appendOutputQname(bb, qname).tab();}
		bb.append(flag).tab();
		if(RNAME_AS_BYTES){
			assert(!(rname==null && rnameS!=null));
			appendOutputRname(bb, rname).tab();
		}else{
			assert(!(rname!=null && rnameS==null)) : RNAME_AS_BYTES+", "+rname+", "+rnameS;
			appendOutputRnameS(bb, rnameS).tab();
		}
		bb.append(pos).tab();
		bb.append(mapq).tab();
		if(cigar==null){bb.append('*');}else{bb.append(cigar);}
		bb.tab();
		appendOutputRname(bb, rnext).tab();
		bb.append(pnext).tab();
		bb.append(tlen).tab();
//		int len=bb.length;
		if(mapped() && strand()==Shared.MINUS && FLIP_ON_LOAD){//[SamToBamConverter#001]/flipsam FIX 2026-06-20 (greenlit):
			//+`&& FLIP_ON_LOAD` so SAM-text output mirrors the load-flip. With flipsam=f the seq is already
			//forward-reference (load didn't flip), so we write it as-is (else branch) instead of double-flipping.
			//Symmetric with SamToBamConverter and BamToSamConverter. VALIDATED 2026-06-20 (Furina): flipsam=f and
			//flipsam=t both write spec-correct forward-reference reverse-strand SEQ (1002 phix rev reads, 0 mismatch,
			//byte-identical BAMs, lossless sam->bam->sam round-trip; pre-fix flipsam=f would have stored RC, which it doesn't).
			appendReverseComplemented(bb, seq).tab();
			appendQualReversed(bb, qual);
//			assert(bb.length==len+seq.length+qual.length+1) : bb.length-len;
		}else{
			appendTo(bb, seq).tab();
			appendQual(bb, qual);
//			assert(bb.length==len+seq.length+qual.length+1) : bb.length-len;
		}
		if(optional!=null){
			for(String s : optional){
				bb.tab().append(s);
			}
		}
		return bb;
	}

	/** Returns SAM text using current output policy, with no final newline. */
	@Override
	public String toString(){return toBytes(null).toString();}

	/** Returns a new default-charset encoding of the suffix after QNAME's sixth underscore.
	 * Requires nonnull QNAME; fewer separators return null, and a final separator yields empty bytes.
	 * This is the legacy positional convention, not a CustomHeader parser. */
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

	/** Tests only that CIGAR is nonnull/nonempty and does not start with '*'; no grammar validation. */
	public boolean hasCigar(){return cigar!=null && cigar.length()>0 && cigar.charAt(0)!='*';}

	/** Tests for X or = anywhere in a nonempty CIGAR not starting with '*'; no grammar validation. */
	public boolean hasCigarXE(){
		if(cigar==null || cigar.length()<1 || cigar.charAt(0)=='*'){return false;}
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(c=='=' || c=='X'){return true;}
		}
		return false;
	}

	/** Tests whether read is part of a paired sequencing template.
	 * @return True if FLAG bit 0x1 is set */
	public boolean hasMate(){return (flag&0x1)==0x1;}

	/** Tests whether read pair is properly aligned.
	 * @return True if FLAG bit 0x2 is set */
	public boolean properPair(){return (flag&0x2)==0x2;}

	/** Tests whether this read is mapped.
	 * @return True if read is mapped (FLAG bit 0x4 not set) */
	public boolean mapped(){
		return (flag&0x4)!=0x4;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}

	/** Tests whether mate read is mapped.
	 * @return True if mate is mapped (FLAG bit 0x8 not set) */
	public boolean nextMapped(){
		return (flag&0x8)!=0x8;
//		0x4 fragment unmapped
//		0x8 next fragment in the template unmapped
	}

	/** Returns strand of this read.
	 * @return 0 for plus strand, 1 for minus strand */
	public byte strand(){return ((flag&0x10)==0x10 ? (byte)1 : (byte)0);}

	/** Returns strand of mate read (alias for nextStrand).
	 * @return 0 for plus strand, 1 for minus strand */
	public byte mateStrand(){return nextStrand();}
	/** Returns strand of mate/next read.
	 * @return 0 for plus strand, 1 for minus strand */
	public byte nextStrand(){return ((flag&0x20)==0x20 ? (byte)1 : (byte)0);}

	/** Tests whether this is the first fragment in template.
	 * @return True if FLAG bit 0x40 is set */
	public boolean firstFragment(){return (flag&0x40)==0x40;}

	/** Tests whether this is the last fragment in template.
	 * @return True if FLAG bit 0x80 is set */
	public boolean lastFragment(){return (flag&0x80)==0x80;}

	/** Returns pair number (0 for first fragment, 1 for last).
	 * @return 0 if first fragment, 1 if last fragment, 0 if neither */
	public int pairnum(){return firstFragment() ? 0 : lastFragment() ? 1 : 0;}

	/** Tests whether the secondary-alignment bit (SAM 0x100) is NOT set.
	 * NOTE: this is NOT the SAM "primary" definition, which also excludes supplementary
	 * (0x800); use primary()/nonPrimary() below for the spec-correct test.
	 * @return True if FLAG bit 0x100 is not set */
	public boolean nonSecondary(){return (flag&0x100)==0;}
	/** Tests whether this is the SAM primary (representative) alignment: neither the
	 * secondary (0x100) nor the supplementary (0x800) bit is set.
	 * @return True if neither 0x100 nor 0x800 is set */
	public boolean primary(){return (flag&0x900)==0;}
	/** Tests whether this is a non-primary alignment: secondary (0x100) or supplementary (0x800).
	 * @return True if 0x100 or 0x800 is set */
	public boolean nonPrimary(){return (flag&0x900)!=0;}
	/** Clears the secondary bit when true, sets it when false; leaves supplementary status unchanged.
	 * Therefore true does not guarantee primary() for an already supplementary record.
	 * @param b Requested non-secondary status */
	public void setPrimary(boolean b){
		if(b){
			flag=flag&~0x100;
		}else{
			flag=flag|0x100;
		}
	}
	/** Sets mapped status of this read.
	 * @param b True if mapped, false if unmapped */
	public void setMapped(boolean b){
		if(b){
			flag=flag&~0x4;
		}else{
			flag=flag|0x4;
		}
	}
	/** Sets first fragment flag.
	 * @param b True if this is first fragment in template */
	public void setFirstFragment(boolean b){
		if(b){
			flag=flag|0x40;
		}else{
			flag=flag&~0x40;
		}
	}
	/** Sets strand of this read.
	 * @param strand 0 for plus, 1 for minus */
	public void setStrand(int strand){
		if(strand==1){
			flag=flag|0x10;
		}else{
			assert(strand==0);
			flag=flag&~0x10;
		}
	}

	/** Tests whether read failed quality controls.
	 * @return True if FLAG bit 0x200 is set */
	public boolean discarded(){return (flag&0x200)==0x200;}

	/** Tests whether read is PCR or optical duplicate.
	 * @return True if FLAG bit 0x400 is set */
	public boolean duplicate(){return (flag&0x400)==0x400;}

	/** Tests whether this is a supplementary alignment.
	 * @return True if FLAG bit 0x800 is set */
	public boolean supplementary(){return (flag&0x800)==0x800;}

	/** Returns true when pairedOnSameChrom is false, or TLEN is zero or positive.
	 * Uses pairedOnSameChrom's name-only comparison; does not require the proper-pair bit. */
	public boolean leftmost(){
		if(!pairedOnSameChrom() || tlen==0){return true;}
		return tlen>0;
	}


//	/** Assumes rname is an integer. */
//	public int chrom(){
//		if(Data.GENOME_BUILD<0){return -1;}
//		HashMap sc
//	}

	/** Tests whether read has ambiguous mapping (low MAPQ).
	 * @return True if mapped with MAPQ &lt; 4 */
	public boolean ambiguous(){return mapped() && mapq<4;}

	/** Legacy numeric-chromosome decoder, disabled by an unconditional assertion under -ea.
	 * With assertions disabled it expects byte-name storage and performs unchecked decimal
	 * accumulation after limited endpoint checks; not a general scaffold-name lookup. */
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

	/** Returns pos-1 minus the selected leading H/S counts, without checking mapped status.
	 * The result may be negative for clipping or missing coordinates. */
	public int start(boolean includeSoftClip, boolean includeHardClip){
		int x=countLeadingClip(cigar, includeSoftClip, includeHardClip);//Leading CIGAR clips only; the Read constructor's match-array clip/indel adjustments are separate.
		return pos-1-x;
	}

	/** Returns inclusive start+reference-span-1 for mapped records with usable CIGAR.
	 * Otherwise returns start plus max(sequence-length-1,0), or start when sequence is absent.
	 * A nonnull mapped CIGAR must be nonempty; selected clip contributions come from calcCigarLength. */
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

	/** Uses stop's inclusive endpoint for mapped records with usable CIGAR.
	 * The fallback instead returns start+sequence-length, or -1 for absent sequence.
	 * These branches have different endpoint conventions; nonnull mapped CIGAR must be nonempty. */
	public int stop2(final int start, final boolean includeSoftClip, final boolean includeHardClip){
		if(mapped() && cigar!=null && cigar.charAt(0)!='*'){return stop(start, includeSoftClip, includeHardClip);}
//		return (seq==null ? -1 : start()+seq.length());
		return (seq==null ? -1 : start+seq.length);
	}

	/** Returns numeric identifier for this read.
	 * @return Always returns 0 (placeholder implementation) */
	public long numericId(){return 0;}

	/** Compares active reference names only: equal names or an '=' token on either side return true.
	 * Both null names compare equal. Does not check mate-presence, mapped or proper-pair bits. */
	public boolean pairedOnSameChrom(){
//		assert(false) : (rname==null ? "nullX" : new String(rname))+", "+
//		(rnext==null ? "nullX" : new String(rnext))+", "+Tools.equals(rnext, byteequals)+", "+Arrays.equals(rname, rnext)+"\n"+this;
		if(RNAME_AS_BYTES){
			return Tools.equals(rnext, byteequals) || Arrays.equals(rname, rnext) || (/*pairnum()==1 &&*/ Tools.equals(rname, byteequals));
		}else{
			return Tools.equals(rnext, byteequals) || Tools.equals(rnameS, rnext) || (/*pairnum()==1 &&*/ stringequals.equals(rnameS));
		}
	}

	/** Parses a signed decimal prefix after QNAME's fifth underscore using the legacy layout.
	 * Fewer separators return -1; no digits after an existing separator yields zero.
	 * Stops at a non-digit other than an initial minus; requires nonnull QNAME and does not check overflow. */
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
				if(i==loc+1 && c=='-'){mult=-1;}else{break;}
			}else{
				sum=(sum*10)+(c-'0');
			}
		}
		return sum*mult;
	}


	/** Returns the byte-name length, otherwise string-name length, otherwise one for the missing token. */
	public int rnameLen(){return (rname==null ? rnameS==null ? 1 : rnameS.length() : rname.length);}

	/** Returns the borrowed byte-name reference, possibly null; asserts RNAME_AS_BYTES. */
	public byte[] rname(){
		assert(RNAME_AS_BYTES);
		return rname;
	}
	/** Returns the borrowed mate-reference bytes, possibly null or the shared equals sentinel. */
	public byte[] rnext(){return rnext;}

	/** Canonicalizes a text-field byte name and retains its reference; asserts byte-name mode.
	 * Does not clear a previously stored string representation. */
	public void setRname(byte[] x){
		assert(RNAME_AS_BYTES);
		rname=canonicalize(x);
	}
	/** Canonicalizes mate-reference text bytes, retaining ordinary arrays and sharing the equals sentinel. */
	public void setRnext(byte[] x){rnext=canonicalize(x);}

	/** Canonicalizes a string reference name; asserts string-name mode and does not clear stored bytes. */
	public void setRnameS(String x){
		assert(!RNAME_AS_BYTES);
		rnameS=canonicalize(x);
	}
	/** Canonicalizes mate-reference text, encoding ordinary strings with the default charset.
	 * Null/star becomes null; equals uses shared bytes; an empty string becomes an empty array. */
	public void setRnextS(String x){rnext=canonicalizeB(x);}

	/** Canonicalizes the supplied text token; '*' becomes null. Does not parse or validate CIGAR operations. */
	public void setCigar(String x){cigar=canonicalize(x);}

	/** Canonicalizes sequence text bytes, retaining ordinary arrays; null/empty/star becomes null. */
	public void setSeq(byte[] x){seq=canonicalize(x);}

	/** Retains numeric Phred scores by reference; null or empty arrays mean missing qualities.
	 * Numeric 42 and 61 are scores, not the text-field '*' and '=' sentinels. */
	public void setQual(byte[] x){qual=(x==null || x.length==0 ? null : x);}

	/** Returns the stored string name, otherwise a new US-ASCII decoding of stored bytes, or null. */
	public String rnameS(){return rnameS!=null ? rnameS : rname==null ? null : new String(rname, StandardCharsets.US_ASCII);}
	/** Returns a new US-ASCII decoding of mate-reference bytes, or null. */
	public String rnextS(){return rnext==null ? null : new String(rnext, StandardCharsets.US_ASCII);}

	/** Returns the Character-whitespace prefix of the preferred stored name; requires a nonnull representation. */
	public String rnamePrefix(){return (rnameS!=null ? toPrefix(rnameS) : toPrefix(rname));}


	/** Appends s to the optional list, allocating it if absent; no tag validation or deduplication.
	 * Existing shallow copies share this mutation when they already share the list.
	 * @param s Complete tag string, conventionally XX:T:value */
	public void addOptionalTag(String s){
		if(optional==null){optional=new ArrayList<String>();}
		optional.add(s);
	}

	/**
	 * Finds optional tag with specified prefix.
	 * @param prefix Tag prefix to search for (e.g., "MD:Z:")
	 * @return First matching tag string or null if not found
	 */
	public String findTag(String prefix){
		if(optional==null){return null;}
		for(String s : optional){
			if(s.startsWith(prefix)){return s;}
		}
		return null;
	}

	/** Returns MD tag string if present.
	 * @return MD tag or null if not present */
	public String mdTag(){return findTag("MD:Z:");}
	/** Returns YQ (mate quality) tag string if present.
	 * @return YQ tag or null if not present */
	public String mateqTag(){return findTag("YQ:i:");}
	/**
	 * Parses the first matching optional tag beginning at fixed character offset five.
	 * An optional leading plus is skipped; otherwise the existing signed-decimal parser
	 * stops at the first nondigit. The offset does not depend on prefix length.
	 * @param prefix Prefix for findTag, normally a complete XX:T: prefix
	 * @return Integer value or Integer.MIN_VALUE if not found
	 */
	public int parseIntFlag(String prefix){
		String tag=findTag(prefix);
		if(tag==null){return Integer.MIN_VALUE;}
		//STR386: SAM integer values may have a leading plus.
		return Parse.parseInt(tag, tag.charAt(5)=='+' ? 6 : 5);
	}
	/**
	 * Parses the complete value of the first matching tag after fixed character offset five.
	 * Uses Java float syntax and rounding, including signs and exponents. The offset
	 * does not depend on prefix length; malformed values are not treated as missing.
	 * @param prefix Prefix for findTag, normally a complete XX:T: prefix
	 * @return Float value or -1 if not found
	 */
	public float parseFloatFlag(String prefix){
		String tag=findTag(prefix);
		//STR386: The prefix-scanning Parse overload truncates valid exponent notation.
		return tag==null ? -1 : Float.parseFloat(tag.substring(5));
	}
	/** Returns mate quality (YQ tag) value.
	 * @return Mate MAPQ value or Integer.MIN_VALUE if not present */
	public int mateq(){return parseIntFlag("YQ:i:");}
	/** Returns mate identity (YJ tag) value.
	 * @return Mate identity percentage or -1 if not present */
	public float mateID(){return parseFloatFlag("YJ:f:");}

	/**
	 * Assigns scafnum from a selected reference name, asserting that scafnum is initially negative.
	 * Uses this record's name when mapped or when a non-sentinel byte RNAME exists; otherwise
	 * tries a non-sentinel RNEXT when the mate-mapped bit is set. With no selected name,
	 * leaves scafnum unchanged. Name lookup retains ScafMap's diagnostics.
	 * @param scafMap Map required when a name is selected
	 * @return Current scafnum after the optional lookup
	 */
	public int setScafnum(ScafMap scafMap){
		assert(scafnum<0);

		String name=null;
		if(mapped() || (rname!=null && rname!=byteequals && rname!=bytestar)){
			name=rnameS();
		}else if(nextMapped() && rnext!=null && rnext!=byteequals && rnext!=bytestar){
			name=new String(rnext, StandardCharsets.US_ASCII);
		}
		if(name!=null){scafnum=scafMap.getNumber(name);}
		return scafnum;
	}

	/** Returns a partial heap-size estimate for fixed overhead, CIGAR, optional tags and byte names.
	 * Omits sequence, quality, QNAME, raw MD and payload storage; does not measure retained heap. */
	public long countBytes(){
		long sum=76;
		sum+=(cigar==null ? 0 : cigar.length()*2+16);
		sum+=(optional==null ? 0 : optional.size()*32+16);
		sum+=(rname==null ? 0 : rname.length+16);
		sum+=(rnext==null ? 0 : rnext.length+16);
		return sum;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Extracts only the FLAG field from a SAM line byte array.
	 * @param s Byte array containing SAM line
	 * @return FLAG value, or -1 if header line
	 */
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

	/**
	 * Extracts only the QNAME field from a SAM line byte array.
	 * @param s Byte array containing SAM line
	 * @return QNAME string, or null if header line or missing
	 */
	public static final String parseNameOnly(byte[] s){
		assert(s!=null && s.length>0) : "Blank line.";
		if(s[0]=='@'){return null;}

		int a=0, b=0;

		while(b<s.length && s[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(s);
		String qname=(b==a+1 && s[a]=='*' ? null : new String(s, a, b-a, StandardCharsets.US_ASCII));
		return qname;
	}


	/** Converts an expanded BBTools match array using VERSION-selected 1.3 or 1.4 rules.
	 * Delegates coordinates, clipping, null handling and optional base-length
	 * assertion to the selected converter; does not expand run-length encoded matches. */
	public static String toCigar(byte[] match, int start, int stop, long scafLen, byte[] bases){
		if(SamLine.VERSION>1.3){
			return SamLine.toCigar14(match, start, stop, scafLen, bases);
		}else{
			return SamLine.toCigar13(match, start, stop, scafLen, bases);
		}
	}

	/**
	 * Converts expanded BBTools match operations to CIGAR with M for matches/mismatches.
	 * Maps m/s/S/N/B to M, I/X/Y to I, D to D and C to S. SOFT_CLIP takes precedence
	 * outside the reference, omitting out-of-bounds deletions from soft-clip counts.
	 * Deletion runs longer than INTRON_LIMIT are emitted as N. Input arrays are not modified.
	 * Equal inclusive endpoints describe one reference base and are converted normally.
	 * @param match Expanded match operations; unsupported in-bounds symbols throw
	 * @param readStart Zero-based initial reference position, possibly outside the reference
	 * @param readStop Inclusive stop; retained for API compatibility, not used in conversion
	 * @param reflen Reference length used by the clipping policy
	 * @param bases Optional query array used only for an assertion on the generated query length
	 * @return New CIGAR string, or null for null match
	 */
	public static String toCigar13(byte[] match, int readStart, int readStop, long reflen, byte[] bases){
		if(match==null){return null;}
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
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}else{sb.append(lastMode);}
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
		if(mode=='D' && count>INTRON_LIMIT){sb.append('N');}else{sb.append(mode);}
		if(mode!='D'){cigarlen+=count;}
		opcount+=count;

		assert(bases==null || cigarlen==bases.length) : "\n(cigarlen = "+cigarlen+") != (bases.length = "+(bases==null ? -1 : bases.length)+")\n" +
				"cigar = "+sb+"\nmatch = "+new String(match)+"\nbases = "+new String(bases)+"\n";

		return sb.toString();
	}

	/**
	 * Replaces adjacent M, = and X runs with their combined M run, preserving other operations.
	 * Does not validate operation names or positive lengths; asserts no nonzero trailing count.
	 * @param cigar14 CIGAR text, or null
	 * @return Converted text, null for null, or empty text for empty input
	 */
	public static String toCigar13(String cigar14){
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


	/**
	 * Converts expanded BBTools match operations to CIGAR using =, X and ambiguous M.
	 * Maps m/s to =, S/V to X, N/B to M, I/X/Y to I, D to D and C to S.
	 * SOFT_CLIP takes precedence outside the reference, omitting out-of-bounds deletions
	 * from soft-clip counts. Deletion runs longer than INTRON_LIMIT become N.
	 * Input arrays are not modified; no run-length expansion is performed.
	 * Equal inclusive endpoints describe one reference base and are converted normally.
	 * @param match Expanded match operations; unsupported in-bounds symbols throw
	 * @param readStart Zero-based initial reference position, possibly outside the reference
	 * @param readStop Inclusive stop; retained for API compatibility, not used in conversion
	 * @param reflen Reference length used by the clipping policy
	 * @param bases Optional query array used only for an assertion on the generated query length
	 * @return New CIGAR string, or null for null match
	 */
	public static String toCigar14(byte[] match, int readStart, int readStop, long reflen, byte[] bases){
//		assert(false) : readStart+", "+readStop+", "+reflen;
		if(match==null){return null;}
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
					if(lastMode=='D' && count>INTRON_LIMIT){sb.append('N');}else{sb.append(lastMode);}
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

	/** Sums M/=/X/D/N counts, optionally adding S and H; I contributes zero.
	 * Null returns zero. P and unknown operations throw; this is not a complete CIGAR
	 * validator and an unterminated final digit run is not included in the result. */
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
				}else if(c=='H'){
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

	/** Sums query-consuming M/=/X/I counts, optionally adding S and H; D/N contribute zero.
	 * Null returns zero. P and unknown operations throw; an unterminated final digit run
	 * is not included. Hard clipping contributes only when explicitly requested. */
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
				}else if(c=='H'){
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

	/** Sums query counts with the same operation/clip policy as calcCigarReadLength.
	 * Null returns zero; P/unknown operations throw; a final digit run is not included. */
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
				}else if(c=='H'){
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

	/** Sums selected H/S counts in the initial clipping prefix, stopping at another operation.
	 * Unselected clip types are skipped without ending the prefix. Null or both flags false
	 * returns zero. Expects count/operation text rather than validating its grammar. */
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

	/** Counts selected S/H operations at the end of a valid CIGAR, honoring flags independently.
	 * S may precede a final H. Null, empty, absent clips or both flags false return zero.
	 * Retains the legacy zero result when a selected soft count begins at character zero,
	 * including any accumulated hard count; see the remaining all-clipped question below. */
	public static int countTrailingClip(String cigar, boolean includeSoftClip, boolean includeHardClip){
		if(cigar==null || (!includeSoftClip && !includeHardClip)){return 0;}
		int last=cigar.length()-1;
		final boolean trailingHard=last>=0 && cigar.charAt(last)=='H';
		int len=(includeHardClip && trailingHard ? countTrailingHardClip(cigar) : 0);
		if(!includeSoftClip){return len;}
		if(trailingHard){
			for(last--; last>=0 && Tools.isDigit(cigar.charAt(last)); last--){} //Skip the final H count.
		}
		if(last<0 || cigar.charAt(last)!='S'){return len;}

		int mult=1;
		int i;
		for(i=last-1; i>=0; i--){
			char c=cigar.charAt(i);
			if(Tools.isLetter(c) || c=='='){
				break;
			}
			len+=(c-'0')*mult;
			mult*=10;
		}
		//STR-332 fixes ordinary suffix/flag counts; keep the separate all-clipped policy question.
		//TODO [stream/SamLine#002 remainder]: a soft count starting at zero still returns zero,
		//including accumulated H. calcLeftClip/calcRightClip assert against all-soft clipping;
		//do not change that policy implicitly. countTrailingHardClip also retains its i<0 zero.
		if(i<0){return 0;}
		return len;
	}

	/** Reads the decimal count before the last H, without verifying a trailing suffix.
	 * Null, absent H or a count reaching the string start returns zero; see #002 above. */
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
			len+=(c-'0')*mult;
			mult*=10;
		}
		if(i<0){return 0;}
		return len;
	}

	/**
	 * Counts substitution characters outside deletion runs in an MD payload.
	 * A caret starts deletion mode and a digit ends it; does not validate base symbols.
	 * @param mdTag Nonnull value, with or without the MD:Z: prefix; asserted nonnull
	 * @return Substitution count; null returns zero only when assertions are disabled
	 */
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
						if(mode==NORMAL){normals+=current;}else{assert(false) : mode+", "+current;}
						current=0;
					}
					if(c=='^'){mode=DEL;}else if(mode==DEL){
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

	/** Counts initial C bases in expanded or run-length encoded BBTools match text.
	 * Null/empty or a first operation other than C returns zero; digits encode total run length. */
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

	/** Returns the byte length of the initial C/digit prefix, not its expanded base count.
	 * For example C10m returns three; null/empty or a non-C first byte returns zero. */
	public static int countLeadingClip2(byte[] match){
		if(match==null || match.length<1 || match[0]!='C'){return 0;}
		int mloc=0;
		for(; mloc<match.length; mloc++){
			byte b=match[mloc];
			if(b!='C' && !Tools.isDigit(b)){return mloc;}
		}
		return match.length;
	}

	/** Counts trailing C operations in expanded BBTools match text; null returns zero.
	 * The scanned suffix must not contain run-length digits, checked by assertion. */
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

	/** Returns deletions minus insertions while advancing an out-of-bounds leading match prefix.
	 * Scans expanded operations until reference position reaches zero or match ends. D advances
	 * reference only, I advances query only, and other symbols advance both. Null or rloc>=0
	 * returns zero; scanned run-length digits are rejected by assertion. The result is signed. */
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

	/** Returns deletions minus insertions while retreating through an out-of-bounds trailing suffix.
	 * Scans expanded operations until rloc falls below rlen or match ends. D retreats reference
	 * only, I retreats query only, and other symbols retreat both. Null or rloc below rlen
	 * returns zero; scanned run-length digits are rejected by assertion. The result is signed. */
	public static int countTrailingIndels(int rloc, int rlen, byte[] match){
		//[stream/SamLine#003] FIXED 2026-06-20 (greenlit by Brian): ENABLED the trailing out-of-bounds indel
		//correction (was silently disabled by a copy-pasted leading-version guard). Two coupled changes:
		//(1) guard 'rloc>=0' -> 'rloc<rlen' so it fires on a 3' overhang (rloc>=rlen), mirroring the leading
		//version's rloc<0; (2) loop start 'match.length' -> 'match.length-1' to avoid the AIOOBE (match[length]).
		//Effect: corrects SAM POS/TLEN for reads overhanging a scaffold 3' end with indels in the overhang
		//(previously only the 5' end was corrected). VALIDATED 2026-06-20 (Furina): direct unit test of this method —
		//guard (rloc<rlen->0), del/ins counting (+1/-1), null/empty inputs, AND the deep-overhang case (rloc>>match.length)
		//that the pre-fix loop start 'match.length' would have AIOOBE'd on — all pass; numerically symmetric with the
		//working countLeadingIndels. (Live BBMap 3'-overhang-with-indel repro deferred as impractical per Brian; the
		//unit test exercises the fixed code directly, which is the stronger check.)
		if(match==null || rloc<rlen){return 0;}
		int dels=0;
		int inss=0;
		int cloc=0;
		for(int mloc=match.length-1; mloc>=0 && rloc>=rlen; mloc--){
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

	/** Returns the largest individual operation count in each of five categories.
	 * Slots are M/=, X, D/N, I, and S/H/P. Adjacent operations are not merged;
	 * unknown operations and trailing digits contribute nothing. Null returns null.
	 * @param cigar Count/operation text
	 * @return Newly allocated five-element array, or null; empty input yields all zeros */
	public static final int[] cigarToMdsiMax(String cigar){
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

	/** Sums operation lengths into slots M/=, X, D/N, I, and S/H/P, respectively.
	 * Unknown operations and a trailing digit run contribute nothing; no grammar validation.
	 * @param cigar Count/operation text, or null
	 * @return New five-element totals array, null for null, or all zeros for empty input */
	public static final int[] cigarToMsdic(String cigar){
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

	/** Reclassifies an expanded match array in place using query and reference bases.
	 * In-bounds m/S/N/s (and C when unClip) become m, S or N after case-insensitive comparison.
	 * Out-of-bounds non-indel entries become C; out-of-bounds I/D trigger assertions.
	 * I/X/Y advance query only, D reference only, and aligned/clipped entries both.
	 * Every iteration asserts a valid query cursor, including deletion entries.
	 * @param call Nonnull query bases in the orientation used by match
	 * @param ref Nonnull reference bases
	 * @param match Nonnull expanded operations, modified in place; may be empty
	 * @param refstart Zero-based reference position corresponding to the first match entry
	 * @param unClip Whether to reclassify in-bounds C entries */
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

	/** Converts CIGAR operations directly to short BBTools match text without MD/reference lookup.
	 * Maps =/X/D-or-N/I/S to m/S/D/I/C; M becomes N when allowed, otherwise returns null.
	 * H, P and unknown operations are omitted. Counts follow symbols only when greater than one;
	 * adjacent equal operations are not merged. No complete grammar or length validation.
	 * @param cigar Nonnull CIGAR text; empty text produces an empty array
	 * @param allowM Whether to retain unresolved M as N
	 * @return Newly allocated match bytes, or null upon encountering a disallowed M */
	public static final byte[] cigarToShortMatch_old(String cigar, boolean allowM){

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


	/**
	 * Creates YS:i: with pos plus selected span minus one, without coordinate normalization.
	 * Uses seqLength when perfect or CIGAR is null; otherwise uses reference span including S,
	 * excluding H. The returned coordinate uses the same origin as the supplied pos.
	 * @param pos Start position (normally the SAM POS field)
	 * @param seqLength Query length for the direct path
	 * @param cigar CIGAR for the non-perfect path
	 * @param perfect Select the supplied sequence length
	 * @return Complete YS tag
	 */
	public static String makeStopTag(int pos, int seqLength, String cigar, boolean perfect){
//		return "YS:i:"+(pos+((cigar==null || perfect) ? seqLength : -countLeadingClip(cigar, false)+calcCigarLength(cigar, false))-1); //123456789
		return "YS:i:"+(pos+((cigar==null || perfect) ? seqLength : calcCigarLength(cigar, true, false))-1);
	}

	/**
	 * Creates YL:Z: with two comma-separated lengths. Perfect/null-CIGAR records use seqLength twice.
	 * Otherwise the first value subtracts leading soft clips only from seqLength; the second
	 * is reference span excluding both clip types. The first value retains trailing soft clips.
	 * @param pos Unused
	 * @param seqLength Supplied query length
	 * @param cigar CIGAR for the non-perfect path
	 * @param perfect Select the direct two-length path
	 * @return Complete YL tag
	 */
	public static String makeLengthTag(int pos, int seqLength, String cigar, boolean perfect){
		if(cigar==null || perfect){return "YL:Z:"+seqLength+","+seqLength;}
		return "YL:Z:"+(seqLength-countLeadingClip(cigar, true, false))+","+calcCigarLength(cigar, false, false);
	}

	/**
	 * Creates YI:f:100 when perfect; otherwise formats 100*Read.identity(match) to two decimals.
	 * Uses Read's configured identity policy and Locale.ROOT formatting, not calcIdentity().
	 * @param match Match representation accepted by Read.identity; ignored when perfect
	 * @param perfect Bypass identity calculation with the literal 100 tag
	 * @return Complete identity-percentage tag
	 */
	public static String makeIdentityTag(byte[] match, boolean perfect){
		if(perfect){return "YI:f:100";}
		float f=Read.identity(match);
		return Tools.format("YI:f:%.2f", (100*f));
	}

	/**
	 * Creates YR:i: from the supplied score, without range checking or rescaling.
	 * @param score Alignment score
	 * @return Complete YR tag
	 */
	public static String makeScoreTag(int score){return "YR:i:"+score;}

	/**
	 * Builds MD:Z: from expanded match operations and a loaded Data chromosome.
	 * Uses reference bases for substitutions/deletions and compares raw query/reference bytes
	 * for N operations. C and out-of-scaffold positions are skipped. Does not flip or expand inputs.
	 * Deletions are emitted upon the following non-D operation only when within INTRON_LIMIT;
	 * terminal pending deletions are not flushed. Retains the long-indel limitation below.
	 * Every iteration reads a query base, including deletion entries; callers must provide a
	 * valid query cursor throughout. Input arrays are read, not modified.
	 * @param chrom Loaded chromosome index; a negative index returns null
	 * @param refstart Zero-based chromosome coordinate at the first match entry
	 * @param match Expanded operations, or null to return null
	 * @param call Query bases in the orientation corresponding to match
	 * @param scafloc Inclusive zero-based scaffold start on the chromosome
	 * @param scaflen Scaffold length; upper bound is scafloc+scaflen, exclusive
	 * @return Complete MD tag, or null for null match/negative chromosome
	 */
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
				//STR-326: each deletion run must be tested against INTRON_LIMIT independently.
				//Reset after skipped runs too, or a later short deletion inherits the skipped count
				//and can disappear from MD even though the CIGAR and NM tag retain it.
				if(dels<=INTRON_LIMIT){//Otherwise, ignore it
					md.append(count);
					count=0;
					md.append('^');
					for(int i=rpos-dels; i<rpos; i++){
						md.append((char)cha.get(i));
					}
				}
				dels=0;
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

	/**
	 * Reads the first operation count only, returning it if the operation is S.
	 * Does not skip leading hard clips; asserts that S is not the entire CIGAR.
	 * @param cig CIGAR string
	 * @param id Unused
	 * @return Initial S count, or zero for null CIGAR or another first operation
	 */
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

	/**
	 * Reads the final S count only when S is the last character; does not skip trailing H.
	 * Asserts that a preceding operation exists at an index greater than zero.
	 * @param cig CIGAR string
	 * @param id Read identifier for error reporting
	 * @return Final S count, or zero for null/empty CIGAR or a different final character
	 */
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

//	public int length(boolean includeSoftClip){
//		assert((seq!=null && (seq.length!=1 || seq[0]!='*')) || cigar!=null) :
//			"This program requires bases or a cigar string for every sam line.  Problem line:\n"+this+"\n";
//		return seq==null ? calcCigarBases(cigar, includeSoftClip, false) : seq.length;
//	}

	/**
	 * Uses cached neural MAPQ for a primary Read when ss is null and the selected cache is present.
	 * STRICT_MAPQ selects strict versus loose cache. Otherwise passes r's length/mapping/ambiguity
	 * and the selected alignment score to the numeric overload; no upper clamp is added here.
	 * @param r Nonnull read with alignment data
	 * @param ss Optional site supplying slowScore instead of r.mapScore and bypassing the neural cache
	 * @return Cached nonnegative value or the numeric-overload result
	 */
	public static int toMapq(Read r, SiteScore ss){
		assert(r!=null);
		if(ss==null && r.primary()){
			final int neural=STRICT_MAPQ?NeuralMapqCache.getStrict(r):NeuralMapqCache.getLoose(r);
			if(neural>=0){return neural;}
		}
		int score=(ss==null ? r.mapScore : ss.slowScore);
		return toMapq(score, r.length(), r.mapped(), r.ambiguous());
	}

	/**
	 * Converts alignment score using the legacy length-scaled formula, without an upper cap.
	 * Unmapped/nonpositive-length input returns zero. Penalized ambiguous reads have a floor
	 * of one; the other branch has a floor of four and a logarithmic length factor.
	 * Does not validate score range or guarantee a value at most 255.
	 * @param score Raw alignment score
	 * @param length Query sequence length
	 * @param mapped Whether read is mapped
	 * @param ambig Whether read has ambiguous mapping
	 * @return Formula result with the selected lower bound, or zero for unmapped/empty input
	 */
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

	/** Appends raw bytes, or '*' for null/star placeholders; an empty array emits no bytes.
	 * @param sb Nonnull destination
	 * @param a Borrowed text/base bytes, not numeric qualities
	 * @return The supplied builder */
	private static ByteBuilder appendTo(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		return sb.append(a);
	}

	/** Appends '*' for null; otherwise emits all of a or its prefix under TRIM_QNAME.
	 * The output-prefix helper stops at a character at or below space; empty input stays empty. */
	private static ByteBuilder appendOutputQname(ByteBuilder sb, String a){
		if(a==null){return sb.append('*');}
		if(Shared.TRIM_QNAME){return sb.appendUntilWhitespace(a);}
		return sb.append(a);
	}

	/** Appends '*' for null/star, otherwise raw or TRIM_RNAME-prefixed bytes.
	 * The prefix helper stops at a byte at or below space; empty arrays emit no bytes. */
	private static ByteBuilder appendOutputRname(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}
		if(Shared.TRIM_RNAME){return sb.appendUntilWhitespace(a);}
		return sb.append(a);
	}

	/** Appends '*' for null/star, otherwise all or the TRIM_RNAME prefix of the string.
	 * The prefix helper stops at a character at or below space; empty input stays empty. */
	private static ByteBuilder appendOutputRnameS(ByteBuilder sb, String a){
		if(a==null || a.equals("*")){return sb.append('*');}
		if(Shared.TRIM_RNAME){return sb.appendUntilWhitespace(a);}
		return sb.append(a);
	}

	/** Appends a text string, or '*' for null/star; an empty string emits no characters.
	 * @param sb Nonnull destination
	 * @param a Text to append, without trimming
	 * @return The supplied builder */
	private static ByteBuilder appendTo(ByteBuilder sb, String a){
		if(a==null || a==stringstar || (a.length()==1 && a.charAt(0)=='*')){return sb.append('*');}
		return sb.append(a);
	}

	/** Appends the extended-alphabet reverse complement without modifying a.
	 * Null/star emits '*'; empty arrays emit no bytes. Requires table-indexable input symbols
	 * and independent destination storage.
	 * @param sb Nonnull destination
	 * @param a Borrowed sequence bytes
	 * @return The supplied builder */
	private static ByteBuilder appendReverseComplemented(ByteBuilder sb, byte[] a){
		if(a==null || a==bytestar || (a.length==1 && a[0]=='*')){return sb.append('*');}

		sb.ensureExtra(a.length);
		byte[] buffer=sb.array;
		int i=sb.length;
		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=AminoAcid.baseToComplementExtended[a[j]];}
		sb.length+=a.length;

		return sb;
	}

	/** Appends numeric qualities plus 33 without clamping or modifying the input.
	 * Null/legacy canonical missing arrays emit '*'; empty arrays emit no bytes.
	 * @param sb Nonnull destination with storage independent of a
	 * @param a Borrowed numeric Phred scores
	 * @return The supplied builder */
	private static ByteBuilder appendQual(ByteBuilder sb, byte[] a){
		//Only the legacy missing-array identity is a sentinel; numeric Phred 42 is valid.
		if(a==null || a==bytestar){return sb.append('*');}

//		sb.ensureExtra(a.length);
//		byte[] buffer=sb.array;
//		int i=sb.length;
//		for(int j=0; j<a.length; i++, j++){buffer[i]=(byte)(a[j]+33);}
//		sb.length+=a.length;
		Vector.addAndAppend(a, sb, 33);

		return sb;
	}

	/** Appends numeric qualities in reverse order, adding 33 without clamping.
	 * Does not reverse a in place. Null/legacy canonical missing arrays emit '*'; empty arrays emit nothing.
	 * @param sb Nonnull destination with storage independent of a
	 * @param a Borrowed numeric Phred scores
	 * @return The supplied builder */
	private static ByteBuilder appendQualReversed(ByteBuilder sb, byte[] a){
		//Only the legacy missing-array identity is a sentinel; numeric Phred 42 is valid.
		if(a==null || a==bytestar){return sb.append('*');}

//		sb.ensureExtra(a.length);
//		byte[] buffer=sb.array;
//		int i=sb.length;
//		for(int j=a.length-1; j>=0; i++, j--){buffer[i]=(byte)(a[j]+33);}
//		sb.length+=a.length;
		Vector.addAndAppendReversed(a, sb, 33);

		return sb;
	}



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


	/**
	 * Builds selected SAM bits from read state, without validating an alignment or setting duplicate.
	 * A nonnull r2 sets paired status and enables fragment bits. Proper-pair requires both reads
	 * mapped/valid with match arrays, r.paired(), and sameScaf. Secondary, discarded and supplementary
	 * bits come from r; strand bits are copied independently of mapped status.
	 * @param r Nonnull source read
	 * @param r2 Mate context, possibly null
	 * @param fragNum Zero sets first-fragment; positive sets last-fragment; negative sets neither
	 * @param sameScaf Caller-supplied same-scaffold condition for proper-pair selection
	 * @return SAM FLAG bit field
	 */
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
		if(r.supplementary()){flag|=0x800;}
		return flag;
	}

	/**
	 * Tests whether read is mapped based on FLAG value.
	 * @param flag SAM FLAG bit field
	 * @return True if read is mapped (FLAG bit 0x4 not set)
	 */
	public static boolean mapped(int flag){return (flag&0x4)!=0x4;}

	/**
	 * Extracts strand from FLAG value.
	 * @param flag SAM FLAG bit field
	 * @return 0 for plus strand, 1 for minus strand
	 */
	public static byte strand(int flag){return ((flag&0x10)==0x10 ? (byte)1 : (byte)0);}

	/** Maps '*' to null and '=' to the canonical equals string; preserves other strings, including empty/null. */
	public static final String canonicalize(String x){
		if(x!=null && x.length()>1){return x;}else if(stringstar.equals(x)){return null;}else if(stringequals.equals(x)){return stringequals;}else if(x==null){return x;/* handle? */}else{return x;}
	}

	/** Maps null/star to null and equals to shared byteequals; default-charset encodes other strings.
	 * Empty strings become empty arrays. Treat the canonical equals array as borrowed. */
	public static final byte[] canonicalizeB(String x){
		if(x!=null && x.length()>1){return x.getBytes();}else if(stringstar.equals(x)){return null;}else if(stringequals.equals(x)){return byteequals;}else if(x==null){return null;/* handle? */}else{return x.getBytes();}
	}

	/** Canonicalizes text-field bytes: null/empty/star becomes null and equals uses shared byteequals.
	 * Other arrays are retained. This is not numeric-quality normalization; canonical arrays are borrowed. */
	public static final byte[] canonicalize(byte[] x){
		if(x!=null && x.length>1){return x;}else if(x==null || x.length==0){return null;}else if(x[0]==star){return null;}else if(x[0]==equals){return byteequals;}else{return x;}
	}

	/**
	 * Extracts prefix of string before first whitespace character.
	 * @param s Input string
	 * @return Prefix before whitespace or full string
	 */
	private static String toPrefix(String s){
		for(int i=0; i<s.length(); i++){
			if(Character.isWhitespace(s.charAt(i))){
				return s.substring(0, i);
			}
		}
		return s;
	}

	/**
	 * Extracts prefix of byte array before first whitespace character.
	 * @param s Input byte array
	 * @return Prefix string before whitespace or full string
	 */
	private static String toPrefix(byte[] s){
		for(int i=0; i<s.length; i++){
			if(Character.isWhitespace(s[i])){
				return new String(s, 0, i, StandardCharsets.US_ASCII);
			}
		}
		return new String(s, StandardCharsets.US_ASCII);
	}


	/** Tests whether any read-group metadata field or READGROUP_TAG is nonnull.
	 * Ignores NO_TAGS and does not ensure READGROUP_ID/READGROUP_TAG are mutually consistent. */
	public static boolean makeReadgroupTags(){
		return READGROUP_ID!=null || READGROUP_CN!=null || READGROUP_DS!=null || READGROUP_DT!=null ||
				READGROUP_FO!=null || READGROUP_KS!=null || READGROUP_LB!=null || READGROUP_PG!=null ||
				READGROUP_PI!=null || READGROUP_PL!=null || READGROUP_PU!=null || READGROUP_SM!=null ||
				READGROUP_TAG!=null;
	}

	/** Returns a historical subset-of-flags summary, gated by NO_TAGS.
	 * Omits MAKE_MD_TAG/MAKE_XT_TAG and includes MAKE_AS_TAG, which makeOptionalTags does not use.
	 * This is not a prediction of whether a particular record will receive optional tags. */
	public static boolean makeOtherTags(){
		if(NO_TAGS){return false;}
		return MAKE_AM_TAG || MAKE_NM_TAG || MAKE_SM_TAG || MAKE_XM_TAG || MAKE_XS_TAG || MAKE_AS_TAG ||
				MAKE_NH_TAG || MAKE_TOPHAT_TAGS || MAKE_IDENTITY_TAG || MAKE_SCORE_TAG || MAKE_STOP_TAG || MAKE_LENGTH_TAG ||
				MAKE_CUSTOM_TAGS || MAKE_INSERT_TAG || MAKE_CORRECTNESS_TAG || MAKE_TIME_TAG || MAKE_BOUNDS_TAG || MAKE_MATEQ_TAG || MAKE_DUAL_MAPQ_TAGS;
	}

	/** ORs the read-group metadata and historical other-tag summaries.
	 * Read-group metadata can make this true even with NO_TAGS; not a generation guarantee. */
	public static boolean makeAnyTags(){return makeReadgroupTags() || makeOtherTags();}

	/*--------------------------------------------------------------*/
	/*--------------------        Fields        --------------------*/
	/*--------------------------------------------------------------*/

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



	/** Serialization version identifier */
	private static final long serialVersionUID=-4180486051387471116L;

	/** Query template name (QNAME field) */
	public String qname;
	/** Bitwise FLAG field containing read pair and mapping information */
	public int flag;
	/** 1-based leftmost mapping position (POS field) */
	public int pos;
	/** Mapping quality score (MAPQ field) */
	public int mapq;
	/** CIGAR string describing alignment operations */
	public String cigar;
	/** Position of mate/next read (PNEXT field) */
	public int pnext;
	/** Observed template length (TLEN field) */
	public int tlen;
	/** Mutable sequence bytes; orientation follows construction/load policy and may be shared. */
	public byte[] seq;
	/** Numeric Phred quality scores, or null when absent; SAM text output adds 33. */
	public byte[] qual;
	/** Mutable optional-tag list, shared by shallow copies; entries retain their full tag prefixes. */
	public ArrayList<String> optional;
	/** Optional raw MD payload used by callers; mdTag() itself searches the optional-tag list. */
	public byte[] mdTag;

	/** Caller-specific payload retained by shallow copies. */
	public Object obj;
	/** Scaffold number for coordinate mapping */
	public int scafnum=-1;


	/** Reference sequence name as byte array when RNAME_AS_BYTES=true */
	private byte[] rname;
	/** Reference name of mate/next read as byte array */
	private byte[] rnext;

	/** Reference sequence name as String when RNAME_AS_BYTES=false */
	private String rnameS;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields        -----------------*/
	/*--------------------------------------------------------------*/

	/** Text-field missing-token byte, not a numeric quality sentinel by value. */
	private static final byte star=(byte)'*';
	/** Text-field same-reference token byte. */
	private static final byte equals=(byte)'=';

	/** Constant for missing string fields in SAM format */
	private static final String stringstar="*";
	/** Constant indicating same reference as mate */
	private static final String stringequals="=";
	/** Shared legacy missing-token array; callers of numeric-quality helpers use identity only. */
	private static final byte[] bytestar=new byte[]{star};
	/** Shared mutable array representing '='; borrowed by records and never to be modified by callers. */
	public static final byte[] byteequals=new byte[]{equals};
	/** Shared XS strand-tag strings selected by the configured library convention. */
	private static final String XSPLUS="XS:A:+", XSMINUS="XS:A:-";
//	private static final double inv100=0.01d;
//	private static float minratio=0.4f;

	/** One-shot chrom_old warning gate, initially true when user.dir contains /bushnell/. */
	private static boolean warning=System.getProperty("user.dir").contains("/bushnell/");

	/** Read group identifier */
	public static String READGROUP_ID=null;
	/** Read group sequencing center name */
	public static String READGROUP_CN=null;
	/** Read group description */
	public static String READGROUP_DS=null;
	/** Read group date */
	public static String READGROUP_DT=null;
	/** Read group flow order */
	public static String READGROUP_FO=null;
	/** Read group key sequence */
	public static String READGROUP_KS=null;
	/** Read group library */
	public static String READGROUP_LB=null;
	/** Read group programs used for processing */
	public static String READGROUP_PG=null;
	/** Read group predicted median insert size */
	public static String READGROUP_PI=null;
	/** Read group platform/technology */
	public static String READGROUP_PL=null;
	/** Read group platform unit */
	public static String READGROUP_PU=null;
	/** Read group sample */
	public static String READGROUP_SM=null;

	/** Complete RG tag appended by makeOptionalTags when READGROUP_ID is nonnull; configured externally. */
	public static String READGROUP_TAG=null;

	/** Enable MD generation from loaded Data chromosome/match arrays. Turn this off for RNAseq or long indels. */
	public static boolean MAKE_MD_TAG=false;

	/** Suppress makeOptionalTags entirely; does not remove stored tags or prevent their serialization. */
	public static boolean NO_TAGS=false;

	/** Generate AM from current MAPQ and the mate-score/length calculation in makeOptionalTags. */
	public static boolean MAKE_AM_TAG=true;
	/** Generate NM from expanded match operations and clip/intron policies, or zero when perfect. */
	public static boolean MAKE_NM_TAG=true;
	/** Copy this record's MAPQ into SM for mapped reads. */
	public static boolean MAKE_SM_TAG=false;
	/** Generate the historical XM site-count estimate unless MAKE_TOPHAT_TAGS selects its own branch. */
	public static boolean MAKE_XM_TAG=false;
	/** Enable XS for mapped records whose CIGAR contains N, using the configured library convention. */
	public static boolean MAKE_XS_TAG=false;
	/** Emit XT:A:R for mapped, ambiguous, non-secondary reads; no unique-read XT is generated here. */
	public static boolean MAKE_XT_TAG=true;
	/** Reserved AS option included in makeOtherTags; makeOptionalTags does not consult this flag. */
	public static boolean MAKE_AS_TAG=false; //TODO: Alignment score from aligner
	/** Emit NH from site-list size when secondary output is enabled with multiple sites, otherwise one. */
	public static boolean MAKE_NH_TAG=false;
	/** Emit historical TopHat compatibility placeholders, taking precedence over the separate XM branch. */
	public static boolean MAKE_TOPHAT_TAGS=false;
	/** Invert the pair-adjusted XS strand selected by makeXSTag. */
	public static boolean XS_SECONDSTRAND=false;
	/** Generate YI from Read.identity (or literal 100 when perfect), under the builder's data-availability guards. */
	public static boolean MAKE_IDENTITY_TAG=false;
	/** Generate YR (raw alignment score) tags */
	public static boolean MAKE_SCORE_TAG=false;
	/** Generate YS (stop position) tags */
	public static boolean MAKE_STOP_TAG=false;
	/** Generate YL lengths using makeLengthTag's leading-soft-clip and reference-span conventions. */
	public static boolean MAKE_LENGTH_TAG=false;
	/** Generate BBTools custom tags (X1, X2, X3, X5, X6, X7) */
	public static boolean MAKE_CUSTOM_TAGS=false;
	/** Generate X8 (insert size) tags */
	public static boolean MAKE_INSERT_TAG=false;
	/** Generate X9 (correctness) tags */
	public static boolean MAKE_CORRECTNESS_TAG=false;
	/** Generate X0 from r.obj, which makeOptionalTags asserts is a nonnull Long. */
	public static boolean MAKE_TIME_TAG=false;
	/** Generate XB from caller-supplied bounds flags; does not independently bypass the unmapped early return. */
	public static boolean MAKE_BOUNDS_TAG=false;
	/** Generate YQ/YJ (mate quality/identity) tags */
	public static boolean MAKE_MATEQ_TAG=false;
	/** Emit cached QL/QS for mapped primary reads; asserts both caches are present or both absent. */
	public static boolean MAKE_DUAL_MAPQ_TAGS=false;
	/** Select strict instead of loose neural MAPQ when toMapq(Read, null) can use a primary read's cache. */
	public static boolean STRICT_MAPQ=false;

	/** Reduce MAPQ for ambiguously mapping reads */
	public static boolean PENALIZE_AMBIG=true;
	/** Enable missing-match conversion in toRead; a CIGAR containing '=' also triggers conversion when false. */
	public static boolean CONVERT_CIGAR_TO_MATCH=true;
	/** Convert out-of-reference match operations to soft clips in toCigar13/14; omit out-of-bounds deletions. */
	public static boolean SOFT_CLIP=true;
	/** Omit shared SEQ/QUAL arrays when rebuilding a secondary Read as SamLine; serialization emits '*'. */
	public static boolean SECONDARY_ALIGNMENT_ASTERISKS=true;
	/** OK to use the "setFrom" function which uses the old SamLine instead of translating the read, if a genome is not loaded. */
	public static boolean SET_FROM_OK=false;
	/** Preserve paired names in cooperating writers and retain simple mate suffixes in the Read constructor. */
	public static boolean KEEP_NAMES=false;
	/** Select CIGAR generation policy: values above 1.3 use =/X-aware toCigar14, otherwise toCigar13. */
	public static float VERSION=1.4f;
	/** CIGAR deletion runs strictly longer than this become N; also consulted by MD/NM and identity helpers. */
	public static int INTRON_LIMIT=Integer.MAX_VALUE;
	/** Select active RNAME byte/string storage; keep consistent with record construction and output policy. */
	public static boolean RNAME_AS_BYTES=true;//Effect on speed is negligible for pileup...

	/** Prefer available MD text over reference bases during toShortMatch correction. */
	public static boolean PREFER_MDTAG=false;
	/** Enable optional X-to-no-call correction when CIGAR contains X but no M.
	 * May require reference lookup and expanded match correction in toShortMatch. */
	public static boolean FIX_MATCH_NS=false;

	/** Records that Parser saw an explicit XS option, even one disabling XS; used by callers to choose defaults. */
	public static boolean setxs=false;
	/** Records an explicit intron-length option, allowing callers to preserve the requested threshold. */
	public static boolean setintron=false;

	/** Request sorted scaffold headers in SamHeader's supporting paths, historically for TopHat compatibility. */
	public static boolean SORT_SCAFFOLDS=false;

	/** Parse QNAME field 0 in the LineParser constructor; false leaves it null. */
	public static boolean PARSE_0=true;
	/** Parse RNAME field 2 into the selected byte/string representation. */
	public static boolean PARSE_2=true;
	/** Parse CIGAR field 5; false leaves it null. */
	public static boolean PARSE_5=true;
	/** Parse RNEXT field 6 into bytes; false leaves it null. */
	public static boolean PARSE_6=true;
	/** Parse PNEXT field 7; false leaves zero. */
	public static boolean PARSE_7=true;
	/** Parse TLEN field 8; false leaves zero. */
	public static boolean PARSE_8=true;
	/** Parse QUAL field 10 and convert from ASCII-33; false leaves numeric qualities null. */
	public static boolean PARSE_10=true;
	/** Parse optional fields after the eleven mandatory columns, subject to the prefix-selection flags. */
	public static boolean PARSE_OPTIONAL=true;
	/** Select MD:-prefixed optional fields; takes precedence over the YQ-only flag. */
	public static boolean PARSE_OPTIONAL_MD_ONLY=false;
	/** Select YQ:-prefixed optional fields when optional parsing is enabled and MD-only selection is false. */
	public static boolean PARSE_OPTIONAL_MATEQ_ONLY=false;
	/** Permit ambiguous M-only cigars (M with no =/X, e.g. minimap2 output) in identity/sub calculations.
	 * Default false: calcIdentity()/countSubs() require a usable MD/reference conversion for such cigars,
	 * rather than guessing their match/mismatch split from CIGAR alone. True permits the fallback of
	 * treating M as matches/zero substitutions when no match array is obtained. A tool that only needs
	 * coverage and not identity (QuickBin) sets this true to process M reads permissively. Note: M ALONGSIDE
	 * =/X is always fine (BBMap marks N-base regions with M) and never triggers the crash. */
	public static boolean M_CIGARS_OK=false;

	/** Flip mapped minus-strand sequence/quality on parsing and reverse that policy for text output.
	 * Configure consistently across record lifetime; toShortMatch's temporary flip is independent. */
	public static boolean FLIP_ON_LOAD=true;

	/** Enable verbose debug output */
	public static boolean verbose=false;

}
