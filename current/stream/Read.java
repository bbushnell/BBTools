package stream;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import align2.GapTools;
import align2.QualityTools;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import hiseq.IlluminaHeaderParser2;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import simd.Vector;
import structures.ByteBuilder;
import structures.FloatList;
import ukmer.Kmer;

/**
 * Mutable nucleotide or amino-acid sequence, qualities, alignment state and caller metadata.
 * Constructors retain supplied arrays; configured validation may modify them immediately.
 * A mate link, the paired-alignment flag and the read1/read2 bit are separate state.
 * Most mutators update only their named fields, not every dependent alignment or flag.
 * Callers own array/object sharing and synchronization. Static validation and formatting
 * policy is process-wide; configure it before concurrent use.
 *
 * @author Brian Bushnell
 */
public final class Read implements Comparable<Read>, Cloneable, Serializable{

	/*--------------------------------------------------------------*/
	/*---------------------        Main        ---------------------*/
	/*--------------------------------------------------------------*/

	/** Prints the supplied match string, its short encoding and two successive long expansions. */
	public static void main(String[] args){
		byte[] a=args[0].getBytes();
		System.out.println(new String(a));
		byte[] b=toShortMatchString(a);
		System.out.println(new String(b));
		byte[] c=toLongMatchString(b);
		System.out.println(new String(c));
		byte[] d=toLongMatchString(c);
		System.out.println(new String(d));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a read with bases, quality scores, and numeric ID.
	 * @param bases_ Sequence bases
	 * @param quals_ Quality scores (may be null)
	 * @param id_ Numeric identifier
	 */
	public Read(byte[] bases_, byte[] quals_, long id_){this(bases_, quals_, Long.toString(id_), id_);}

	/**
	 * Creates a read with bases, quality scores, name, and numeric ID.
	 *
	 * @param bases_ Sequence bases
	 * @param quals_ Quality scores (may be null)
	 * @param name_ Read identifier string
	 * @param id_ Numeric identifier
	 */
	public Read(byte[] bases_, byte[] quals_, String name_, long id_){this(bases_, quals_, name_, id_, VALIDATE_IN_CONSTRUCTOR);}

	/** Retains bases/qualities and optionally validates them in place.
	 * Does not copy arrays or infer mate/alignment state. */
	public Read(byte[] bases_, byte[] quals_, String name_, long id_, boolean validate){
		bases=bases_;
		quality=quals_;
		id=name_;
		numericID=id_;
		
		if(validate){validate(true);}
	}

	/** Retains arrays and the supplied flags, optionally invoking validate(true).
	 * Unlike the coordinate constructor, this overload does not clear VALIDATEDMASK first. */
	public Read(byte[] bases_, byte[] quals_, String name_, long id_, int flag_, boolean validate){
		bases=bases_;
		quality=quals_;
		id=name_;
		numericID=id_;
		flags=flag_;
		
		if(validate){validate(true);}
	}

	/** Retains arrays/flags with unknown coordinates and configured constructor validation. */
	public Read(byte[] bases_, byte[] quals_, String name_, long id_, int flag_){this(bases_, quals_, name_, id_, flag_, -1, -1, -1);}

	/** Retains arrays and alignment coordinates, using the numeric ID as the name. */
	public Read(byte[] s_, byte[] quals_, long id_, int chrom_, int start_, int stop_, byte strand_){this(s_, quals_, Long.toString(id_), id_, (int)strand_, chrom_, start_, stop_);}

	/** Retains arrays and coordinates, clears the supplied validation bit, then optionally validates.
	 * The strand occupies the low flag bit, so a strand code alone is also a valid flags_ value.
	 * Coordinates and flags are stored independently; this does not infer mapping success. */
	public Read(byte[] bases_, byte[] quals_, String id_, long numericID_, int flags_, int chrom_, int start_, int stop_){
		flags=flags_&~VALIDATEDMASK;
		bases=bases_;
		quality=quals_;
		id=id_;
		numericID=numericID_;
		
		chrom=chrom_;
		start=start_;
		stop=stop_;
		
		if(VALIDATE_IN_CONSTRUCTOR){validate(true);}
	}

	/*--------------------------------------------------------------*/
	/*-------------------        Methods        --------------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Applies configured sequence/quality/header checks and records that validation ran.
	 * May mutate arrays, replace/drop qualities, and set junk/discarded flags. A null
	 * base array clears qualities. The validated bit is not a promise that junk is absent.
	 * The helper's passesJunk result is not returned; inspect flags for retained bad input.
	 *
	 * @param processAssertions Enables explicit fatal diagnostics in several helpers;
	 * ordinary Java assertions still depend on JVM assertion enablement
	 * @return True on normal return, including flag-only junk handling; failures may terminate
	 */
	public boolean validate(final boolean processAssertions){
		assert(!validated());
//		assert(false);
		//		if(false){//This causes problems with error-corrected PacBio reads.
		//			boolean x=(quality==null || quality.length<1 || quality[0]<=80 || !FASTQ.DETECT_QUALITY || FASTQ.IGNORE_BAD_QUALITY);
		//			if(!x){
		//				if(processAssertions){
		//					KillSwitch.kill("Quality value ("+quality[0]+") appears too high.\n"+Arrays.toString(quality)+
		//							"\n"+Arrays.toString(bases)+"\n"+numericID+"\n"+id+"\n"+FASTQ.ASCII_OFFSET);
		//				}
		//				return false;
		//			}
		//		}

		if(bases==null){
			quality=null; //I could require this to be true
			if(FIX_HEADER){fixHeader(processAssertions);}
			setValidated(true);
			return true;
		}

		validateQualityLength(processAssertions);

		final boolean passesJunk;

		//		assert(false) : SKIP_SLOW_VALIDATION+","+VALIDATE_BRANCHLESS+","+JUNK_MODE;
		if(SKIP_SLOW_VALIDATION){
			passesJunk=true;
		}else if(!aminoacid()){//Nucleotide path: SIMD (if available) -> branchless -> scalar. NOTE: VALIDATE_BRANCHLESS is `final true`, so the scalar validateCommonCase() branch is DEAD; the live paths are validateCommonCaseVec (SIMD) and validateCommonCase_branchless.
			if(U_TO_T){uToT();}
			if(VALIDATE_VECTOR && Shared.SIMD){
				passesJunk=validateCommonCaseVec(processAssertions);
			}else if(VALIDATE_BRANCHLESS){
				passesJunk=validateCommonCase_branchless(processAssertions);
			}else{
				passesJunk=validateCommonCase(processAssertions);
			}
		}else{
			//			if(U_TO_T){uToT();}  This is amino, so...
			fixCase();
			passesJunk=validateJunk(processAssertions);
			if(CHANGE_QUALITY){fixQuality();}
		}

		if(FIX_HEADER){fixHeader(processAssertions);}

		setValidated(true);

		return true;
	}

	/** Converts U/u in bases in place through Vector's scalar/SIMD dispatch; null is a no-op. */
	private void uToT(){Vector.uToT(bases);}

	/** Applies the T-to-U lookup in place; null bases are unchanged, other bytes must index the table. */
	private void tToU(){
		if(bases==null){return;}
		for(int i=0; i<bases.length; i++){
			bases[i]=AminoAcid.tToU[bases[i]];
		}
	}

	/** Checks extended symbols using the read's nucleotide or amino alphabet and JUNK_MODE.
	 * IGNORE returns true immediately; FIX replaces rejected bytes with N/X without touching qualities.
	 * FLAG marks junk and returns false on the first rejection; other modes may terminate when
	 * processAssertions is true, otherwise mark junk and return false. Requires nonnull table-indexable bases. */
	private boolean validateJunk(boolean processAssertions){
		assert(bases!=null);
		if(JUNK_MODE==IGNORE_JUNK){return true;}
		final byte nocall;
		final byte[] toNum;
		final boolean aa=aminoacid();
		if(aa){
			nocall='X';
			toNum=AminoAcid.acidToNumberExtended;
		}else{
			nocall='N';
			toNum=AminoAcid.baseToNumberExtended;
		}
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int num=toNum[b];
			//			System.err.println(Character.toString(b)+" -> "+num);
			if(num<0){
				if(JUNK_MODE==FIX_JUNK){
					bases[i]=nocall;
				}else if(JUNK_MODE==FLAG_JUNK){
					setJunk(true);
					return false;
				}else{
					if(processAssertions){
						KillSwitch.kill("\nAn input file appears to be misformatted:\n"
							+ "The character with ASCII code "+b+" appeared where "+(aa ? "an amino acid" : "a base")+" was expected"
							+ (b>31 && b<127 ? ": '"+Character.toString((char)b)+"'\n" : ".\n")
							+ "Sequence #"+numericID+"\n"
							+ "Sequence ID '"+id+"'\n"
							+ "Sequence: "+Tools.toStringSafe(bases)+"\n"
							+ "Flags: "+Long.toBinaryString(flags)+"\n\n"
							+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
					}
					setJunk(true);
					return false;
				}
			}
		}
		return true;
	}

	/** Handles unequal base/quality lengths with TOSS, NULLIFY, then FLAG option precedence.
	 * TOSS clears qualities and marks discarded/junk; NULLIFY clears qualities and marks junk;
	 * FLAG retains the mismatched array and marks junk. With none selected, the fatal diagnostic
	 * requires both processAssertions and enabled Java assertions. Assumes nonnull bases. */
	private void validateQualityLength(boolean processAssertions){
		if(quality==null || quality.length==bases.length){return;}
		if(TOSS_BROKEN_QUALITY){
			quality=null;
			setDiscarded(true);
			setJunk(true);
		}else if(NULLIFY_BROKEN_QUALITY){
			quality=null;
			setJunk(true);
		}else if(FLAG_BROKEN_QUALITY){
			setJunk(true);
		}else{//No broken-quality flag set: crash loud on the bases/quals length mismatch (malformed input).
			boolean x=false;
			assert(x=processAssertions);//Brian's runtime-assertion-detection idiom: x becomes processAssertions ONLY under -ea (the assert is evaluated). So the kill below fires under -ea (BBTools default) and is skipped under -da. Cf. Shared.EA().
			if(x){
				KillSwitch.kill("\nMismatch between length of bases and qualities for read "+numericID+" (id="+id+").\n"+
					"# qualities="+quality.length+", # bases="+bases.length+"\n\n"+
					FASTQ.qualToString(quality)+"\n"+new String(bases)+"\n\n"
					+ "This can be bypassed with the flag 'tossbrokenreads' or 'nullifybrokenquality'");
			}
		}
	}

	/** Caps called-symbol qualities and zeros ambiguous symbols using this read's alphabet.
	 * Uses acidToNumber for amino reads, including its stop-symbol policy; nucleotide
	 * definedness would silently erase valid protein qualities (STR-018).
	 * Does nothing for absent qualities or CHANGE_QUALITY=false. Callers must supply
	 * nonnegative Phred scores and a base for every quality. Package access lets FASTQ
	 * normalize amino reads even when constructor validation is disabled.
	 */
	final void fixQuality(){
		if(quality==null || !CHANGE_QUALITY){return;}
		assert(bases!=null && bases.length>=quality.length) :
			"Alphabet-aware quality normalization needs a symbol for every score; read="+id;
		final byte[] toNumber=aminoacid() ? AminoAcid.acidToNumber : AminoAcid.baseToNumber;
		for(int i=0; i<quality.length; i++){
			final byte b=bases[i], q=quality[i];
			quality[i]=(b>=0 && toNumber[b]>=0 ? qMap[q] : 0);
		}
	}

	/** Applies configured case and dot/dash/X maps in place using this read's alphabet.
	 * Uppercasing takes precedence over lowercase-to-nocall; symbol remapping precedes case mapping.
	 * Null bases or no selected conversion are no-ops; bytes must index the selected tables. */
	private void fixCase(){
		if(bases==null || (!DOT_DASH_X_TO_N && !TO_UPPER_CASE && !LOWER_CASE_TO_N)){return;}
		final boolean aa=aminoacid();

		final byte[] caseMap, ddxMap;
		if(!aa){
			caseMap=TO_UPPER_CASE ? AminoAcid.toUpperCase : 
				LOWER_CASE_TO_N ? AminoAcid.lowerCaseToNocall : null;
			ddxMap=DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocall : null;
		}else{
			caseMap=TO_UPPER_CASE ? AminoAcid.toUpperCase : 
				LOWER_CASE_TO_N ? AminoAcid.lowerCaseToNocallAA : null;
			ddxMap=DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocallAA : null;
		}

		//		assert(false) : (AminoAcid.toUpperCase==caseMap)+", "+ddxMap;

		if(DOT_DASH_X_TO_N){
			if(TO_UPPER_CASE || LOWER_CASE_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=caseMap[ddxMap[b]];
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=ddxMap[b];
				}
			}
		}else{
			if(TO_UPPER_CASE || LOWER_CASE_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					bases[i]=caseMap[b];
				}
			}else{
				assert(false);
			}
		}
	}

	/** Applies nucleotide Vector conversions, with lowercase-to-N taking precedence over uppercase.
	 * Returns Vector.toUpperCase's status only for that branch, otherwise true, including no work.
	 * This return is not an independent junk classification. */
	private boolean fixCaseVec(){
		if(bases==null || (!DOT_DASH_X_TO_N && !TO_UPPER_CASE && !LOWER_CASE_TO_N)){return true;}
		if(DOT_DASH_X_TO_N){Vector.dotDashXToN(bases);}
		if(LOWER_CASE_TO_N){Vector.lowerCaseToN(bases);}else if(TO_UPPER_CASE){return Vector.toUpperCase(bases);}
		return true;
	}

	/** Applies vector-dispatched nucleotide normalization and configured junk handling in place.
	 * Normal return is true even when FLAG marks junk. CRASH delegates to the branchless helper
	 * with diagnostics enabled, irrespective of this method's processAssertions argument.
	 * A failed case conversion is treated as a problem before nucleotide classification. */
	private boolean validateCommonCaseVec(boolean processAssertions){
		assert(!aminoacid());
		assert(bases!=null);
		
		boolean problem=false;
		boolean needCap=false;
		if(IUPAC_TO_N || JUNK_MODE==FIX_JUNK_AND_IUPAC){Vector.iupacToN(bases); needCap=true;}
		//Case conversion reports success; negate it before combining with nucleotide validity.
		//Otherwise FLAG_JUNK marks clean converted reads as junk.
		if(TO_UPPER_CASE || LOWER_CASE_TO_N || DOT_DASH_X_TO_N){problem=!fixCaseVec(); needCap=true;}
		if(quality!=null && CHANGE_QUALITY && needCap){
			Vector.capQuality(quality, bases);
			needCap=false;
		}
		
		if(JUNK_MODE==IGNORE_JUNK){return true;}
		problem=problem||!Vector.isNucleotide(bases);
		if(!problem){return true;}
		if(JUNK_MODE==FLAG_JUNK){setJunk(true);}else if(JUNK_MODE==FIX_JUNK){
			Vector.iupacToN(bases);
			Vector.capQuality(quality, bases);
		}else if(JUNK_MODE==CRASH_JUNK){return validateCommonCase_branchless(true);}
		return true;
	}

	/** Normalizes configured nucleotide symbols/qualities, accumulating invalid-symbol sign bits.
	 * Returns true for ignored, accepted or repaired input; flags junk and returns false when retained.
	 * May terminate for rejected input when diagnostics are enabled. Requires compatible arrays
	 * and table-indexable bytes/qualities. FIX_JUNK_AND_IUPAC rescans only after junk was detected. */
	private boolean validateCommonCase_branchless(boolean processAssertions){

		assert(!aminoacid());
		assert(bases!=null);

		if(TO_UPPER_CASE || LOWER_CASE_TO_N){fixCase();}

		final byte nocall='N';
		final byte[] toNum=AminoAcid.baseToNumber;
		final byte[] map=(DOT_DASH_X_TO_N && IUPAC_TO_N ? AminoAcid.baseToACGTN : 
			DOT_DASH_X_TO_N ? AminoAcid.dotDashXToNocall : IUPAC_TO_N ? AminoAcid.iupacToNocall : null);

		if(JUNK_MODE==IGNORE_JUNK){
			if(quality!=null && CHANGE_QUALITY){
				if(DOT_DASH_X_TO_N || IUPAC_TO_N){
					for(int i=0; i<bases.length; i++){
						final byte b=map[bases[i]];
						final byte q=quality[i];
						final int num=toNum[b];
						bases[i]=b;
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}else{
					for(int i=0; i<bases.length; i++){
						final byte b=bases[i];
						final byte q=quality[i];
						final int num=toNum[b];
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}
			}else if(DOT_DASH_X_TO_N || IUPAC_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=map[bases[i]];
					bases[i]=b;
				}
			}
			return true;
		}

		int junkOr=0;
		//	int iupacOr=0;
		final byte[] toNumE=AminoAcid.baseToNumberExtended;

		//	asdf

		if(DOT_DASH_X_TO_N || IUPAC_TO_N){
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					final byte b=map[bases[i]];
					final byte q=quality[i];
					final int numE=toNumE[b];
					final int num=toNum[b];
					junkOr|=numE;
					//				iupacOr|=num;
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					final byte b=map[bases[i]];
					final int numE=toNumE[b];
					//				final int num=toNum[b];
					junkOr|=numE;
					//				iupacOr|=num;
					bases[i]=b;
				}
			}
		}else{
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final byte q=quality[i];
					final int numE=toNumE[b];
					final int num=toNum[b];
					junkOr|=numE;
					//				iupacOr|=num;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final int numE=toNumE[b];
					//				final int num=toNum[b];
					junkOr|=numE;
					//				iupacOr|=num;
				}
			}
		}

		//	System.err.println(junkOr+", "+JUNK_MODE);

		//STUDIED-PRAISE branchless junk detection: each toNumE[b] is NEGATIVE iff b is an invalid base, so OR-ing them all into junkOr accumulates the sign bit - junkOr>=0 iff EVERY char was valid. The entire clean-read common case exits right here with zero per-char branches; only a dirty read (junkOr<0) pays the re-scan below. [verified: baseToNumberExtended maps junk chars to negatives]
		//Common case
		if(junkOr>=0){return true;}
		//	if(junkOr>=0 && (JUNK_MODE!=FIX_JUNK_AND_IUPAC || iupacOr>=0)){return true;}
		//	
		//	assert(junkOr<0 || (JUNK_MODE==FIX_JUNK_AND_IUPAC && iupacOr<0));

		//TODO: I could disable VALIDATE_BRANCHLESS here, if it's not final
		//VALIDATE_BRANCHLESS=false;
		if(JUNK_MODE==FIX_JUNK){
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				final int numE=toNumE[b];

				if(numE<0){
					bases[i]=nocall;
					if(quality!=null){quality[i]=0;}
				}
			}
			return true;
		}else if(JUNK_MODE==FLAG_JUNK){
			setJunk(true);
			return false;
		}else if(JUNK_MODE==FIX_JUNK_AND_IUPAC){
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				final byte c=AminoAcid.baseToACGTN[b];
				bases[i]=c;
			}
			return true;
		}else{
			if(processAssertions){
				int i=0;
				for(; i<bases.length; i++){
					if(toNumE[bases[i]]<0){break;}
				}
				byte b=bases[i];
				KillSwitch.kill("\nAn input file appears to be misformatted:\n"
					+ "The character with ASCII code "+b+" appeared where a base was expected"
					+ (b>31 && b<127 ? ": '"+Character.toString((char)b)+"'\n" : ".\n")
					+ "Sequence #"+numericID+"\n"
					+ "Sequence ID: '"+id+"'\n"
					+ "Sequence: '"+Tools.toStringSafe(bases)+"'\n\n"
					+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
			}
			setJunk(true);
			return false;
		}
	}

	//DEAD via validate(): reached only through the `else` after `else if(VALIDATE_BRANCHLESS)`, but VALIDATE_BRANCHLESS is `final true` (the author's ~L443 TODO acknowledges this) - so validateCommonCase_branchless always wins. Kept as the scalar reference.
	//#001-fix [stream/Read#001]: the 4 junk-crash messages below reported bases[1] (hardcoded 2nd base, usually VALID) instead of bases[i] (the actual offending char) - a misleading diagnostic. Latent LOW (dead path; would mislead if VALIDATE_BRANCHLESS were ever un-finalized). Fixed bases[1]->bases[i].
	/** Legacy scalar nucleotide normalization/reference implementation; not selected by validate().
	 * May remap bases and qualities, return false after marking junk, or issue a fatal diagnostic.
	 * Its option handling is historical and is not asserted equivalent to the active implementations. */
	private boolean validateCommonCase(boolean processAssertions){

		assert(!aminoacid());
		assert(bases!=null);

		if(TO_UPPER_CASE || LOWER_CASE_TO_N){fixCase();}

		final byte nocall='N';
		final byte[] toNum=AminoAcid.baseToNumber;
		final byte[] ddxMap=AminoAcid.dotDashXToNocall;

		if(JUNK_MODE==IGNORE_JUNK){
			if(quality!=null && CHANGE_QUALITY){
				if(DOT_DASH_X_TO_N){
					for(int i=0; i<bases.length; i++){
						final byte b=ddxMap[bases[i]];
						final byte q=quality[i];
						final int num=toNum[b];
						bases[i]=b;
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}else{
					for(int i=0; i<bases.length; i++){
						final byte b=bases[i];
						final byte q=quality[i];
						final int num=toNum[b];
						quality[i]=(num>=0 ? qMap[q] : 0);
					}
				}
			}else if(DOT_DASH_X_TO_N){
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					bases[i]=b;
				}
			}
		}else if(DOT_DASH_X_TO_N){
			final byte[] toNumE=AminoAcid.baseToNumberExtended;
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					final byte q=quality[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
									+ "The character with ASCII code "+bases[i]+" appeared where a base was expected.\n"
									+ "Sequence #"+numericID+"\n"
									+ "Sequence ID '"+id+"'\n"
									+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
									+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					final int num=toNum[b];
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=ddxMap[bases[i]];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
									+ "The character with ASCII code "+bases[i]+" appeared where a base was expected.\n"
									+ "Sequence #"+numericID+"\n"
									+ "Sequence ID '"+id+"'\n"
									+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
									+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					bases[i]=b;
				}
			}
		}else{
			final byte[] toNumE=AminoAcid.baseToNumberExtended;
			if(quality!=null && CHANGE_QUALITY){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					final byte q=quality[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							bases[i]=b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
									+ "The character with ASCII code "+bases[i]+" appeared where a base was expected.\n"
									+ "Sequence #"+numericID+"\n"
									+ "Sequence ID '"+id+"'\n"
									+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
									+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}

					final int num=toNum[b];
					bases[i]=b;
					quality[i]=(num>=0 ? qMap[q] : 0);
				}
			}else{
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					final int numE=toNumE[b];

					if(numE<0){
						if(JUNK_MODE==FIX_JUNK){
							bases[i]=b=nocall;
						}else if(JUNK_MODE==FLAG_JUNK){
							setJunk(true);
							return false;
						}else{
							if(processAssertions){
								KillSwitch.kill("\nAn input file appears to be misformatted:\n"
									+ "The character with ASCII code "+bases[i]+" appeared where a base was expected.\n"
									+ "Sequence #"+numericID+"\n"
									+ "Sequence ID '"+id+"'\n"
									+ "Sequence: "+Tools.toStringSafe(bases)+"\n\n"
									+ "This can be bypassed with the flag 'tossjunk', 'fixjunk', or 'ignorejunk'");
							}
							setJunk(true);
							return false;
						}
					}
				}
			}
		}
		return true;
	}

	/** Replaces id with Tools.fixHeader's printable-header result under the configured null-header policy. */
	private final void fixHeader(boolean processAssertions){id=Tools.fixHeader(id, ALLOW_NULL_HEADER, processAssertions);}

	/** Compares equal-length bases under the no-call/mismatch policy of the six-argument overload. */
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch){return isDuplicateByBases(r, nmax, mmax, qmax, banSameQualityMismatch, false);}



	/** Compares the common base prefix, counting each uppercase-N position once against nmax.
	 * Other disagreements count against mmax and fail when both qualities exceed qmax,
	 * or when equal qualities are forbidden. Such positions require both quality arrays.
	 * Equal lengths are asserted even when allowDifferentLength is true; that option only
	 * relaxes the explicit length rejection with assertions disabled. Requires nonnull bases. */
	public boolean isDuplicateByBases(Read r, int nmax, int mmax, byte qmax, boolean banSameQualityMismatch, boolean allowDifferentLength){
		int n=0, m=0;
		assert(r.length()==bases.length) : "Merging different-length reads is supported but seems to be not useful.";
		if(!allowDifferentLength && r.length()!=bases.length){return false;}
		int minLen=Tools.min(bases.length, r.length());
		for(int i=0; i<minLen; i++){
			byte b1=bases[i];
			byte b2=r.bases[i];
			if(b1=='N' || b2=='N'){
				n++;
				if(n>nmax){return false;}
			}else if(b1!=b2){
				m++;
				if(m>mmax){return false;}
				//TODO: Possible bug [stream/Read#005] - quality[i]/r.quality[i] dereferenced with no null guard; a FASTA read has null quality, so dedup-by-bases with mismatches (mmax>0) on FASTA would NPE. LATENT LOW: isDuplicateByBases is UNCALLED tree-wide (verified grep) - dead code. If revived, guard quality!=null.
				if(quality[i]>qmax && r.quality[i]>qmax){return false;}
				if(banSameQualityMismatch && quality[i]==r.quality[i]){return false;}
			}
		}
		return true;
	}

	/**
	 * Tests mapped coordinate/strand agreement; different lengths use the separate helper.
	 * One-end checks use start on plus and stop on minus. Optional alignment comparison
	 * distinguishes encoded D and I/X/Y categories, not exact base matches or decoded runs.
	 * Two perfect flags bypass that comparison; a missing match array also skips it.
	 * Requires nonnull receiver bases and a distinct nonmate read.
	 *
	 * @param r Read to compare against
	 * @param bothEnds Whether to check both start and stop positions
	 * @param checkAlignment Whether to compare alignment strings for similarity
	 * @return true if reads map to same location and are considered duplicates
	 */
	public boolean isDuplicateByMapping(Read r, boolean bothEnds, boolean checkAlignment){
		if(bases.length!=r.length()){
			return isDuplicateByMappingDifferentLength(r, bothEnds, checkAlignment);
		}
		assert(this!=r && mate!=r);
		assert(!bothEnds || bases.length==r.length());
		if(!mapped() || !r.mapped()){return false;}
		//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}

		//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
		////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
		//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
		//		return true;

		if(chrom!=r.chrom || strand()!=r.strand()){return false;}
		if(bothEnds){
			if(start!=r.start || stop!=r.stop){return false;}
		}else{
			if(strand()==Shared.PLUS){
				if(start!=r.start){return false;}
			}else{
				if(stop!=r.stop){return false;}
			}
		}
		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				if(match.length!=r.match.length){return false;}
				for(int i=0; i<match.length; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}

	/** Tests unequal-length reads at the same strand-specific endpoint; bothEnds always returns false.
	 * Requires mapped bits, matching chromosome/strand and a distinct nonmate read.
	 * Optional alignment comparison checks D and I/X/Y categories over the common encoded
	 * prefix only; two perfect flags or either missing match skip that comparison. */
	public boolean isDuplicateByMappingDifferentLength(Read r, boolean bothEnds, boolean checkAlignment){
		assert(this!=r && mate!=r);
		assert(bases.length!=r.length());
		if(bothEnds){return false;}
		//		assert(!bothEnds || bases.length==r.length());
		if(!mapped() || !r.mapped()){return false;}
		//		if(chrom==-1 && start==-1){return false;}
		if(chrom<1 && start<1){return false;}

		//		if(chrom!=r.chrom || strand()!=r.strand() || start!=r.start){return false;}
		////		if(mate==null && stop!=r.stop){return false;} //For unpaired reads, require both ends match
		//		if(stop!=r.stop){return false;} //For unpaired reads, require both ends match
		//		return true;

		if(chrom!=r.chrom || strand()!=r.strand()){return false;}

		if(strand()==Shared.PLUS){
			if(start!=r.start){return false;}
		}else{
			if(stop!=r.stop){return false;}
		}

		if(checkAlignment){
			if(perfect() && r.perfect()){return true;}
			if(match!=null && r.match!=null){
				int minLen=Tools.min(match.length, r.match.length);
				for(int i=0; i<minLen; i++){
					byte a=match[i];
					byte b=r.match[i];
					if(a!=b){
						if((a=='D') != (b=='D')){return false;}
						if((a=='I' || a=='X' || a=='Y') != (b=='I' || b=='X' || b=='Y')){return false;}
					}
				}
			}
		}
		return true;
	}

	/**
	 * Accumulates duplicate counts and optionally updates this read's bases/qualities in place.
	 * Also merges linked mates when this read has a mate, requiring a corresponding r.mate.
	 * mergeVectors implies mergeN. Without receiver qualities, only uppercase-N positions
	 * are filled; with qualities, merging positions requires donor qualities too.
	 * Quality-vector merging boosts agreements (cap 48), favors higher-quality disagreements,
	 * and emits N/Q0 on equal positive qualities; two zero-quality disagreements are retained.
	 * Equal lengths and reciprocal mate/ID relationships are asserted. Match reconciliation
	 * is partial; existing alignment metadata can require realignment after merging.
	 *
	 * @param r Read to merge with
	 * @param mergeVectors Whether to merge quality vectors
	 * @param mergeN Whether to fill N-calls with corresponding bases
	 */
	public void merge(Read r, boolean mergeVectors, boolean mergeN){mergePrivate(r, mergeVectors, mergeN, true);}

	/** Implements duplicate accumulation, optional sequence consensus and one level of mate merging.
	 * Requires nonnull reads/bases, compatible qualities for used positions, and distinct objects.
	 * The historical unequal-length extension is fenced by an assertion; if enabled with -da,
	 * it may grow receiver arrays/coordinates and clear both reads' match arrays. */
	private void mergePrivate(Read r, boolean mergeVectors, boolean mergeN, boolean mergeMate){
		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
		mergeN=(mergeN||mergeVectors);

		assert(r.length()==bases.length) : "Merging different-length reads is supported but seems to be not useful.";

		if((mergeN || mergeVectors) && bases.length<r.length()){
			int oldLenB=bases.length;
			start=Tools.min(start, r.start);
			stop=Tools.max(stop, r.stop);
			mapScore=Tools.max(mapScore, r.mapScore);

			bases=KillSwitch.copyOfRange(bases, 0, r.length());
			quality=KillSwitch.copyOfRange(quality, 0, r.quality.length);
			for(int i=oldLenB; i<bases.length; i++){
				bases[i]='N';
				quality[i]=0;
			}
			match=null;
			r.match=null;
		}

		copies+=r.copies;


		//		if(numericID==11063941 || r.numericID==11063941 || numericID==8715632){
		//			System.err.println("***************");
		//			System.err.println(this.toText()+"\n");
		//			System.err.println(r.toText()+"\n");
		//			System.err.println(mergeVectors+", "+mergeN+", "+mergeMate+"\n");
		//		}

		boolean pflag1=perfect();
		boolean pflag2=r.perfect();

		final int minLenB=Tools.min(bases.length, r.length());

		if(mergeN){
			if(quality==null){
				for(int i=0; i<minLenB; i++){
					byte b=r.bases[i];
					if(bases[i]=='N' && b!='N'){bases[i]=b;}
				}
			}else{
				for(int i=0; i<minLenB; i++){
					final byte b1=bases[i];
					final byte b2=r.bases[i];
					final byte q1=Tools.max((byte)0, quality[i]);
					final byte q2=Tools.max((byte)0, r.quality[i]);
					if(b1==b2){
						if(b1=='N'){
							//do nothing
						}else if(mergeVectors){
							//merge qualities
							//						quality[i]=(byte) Tools.min(40, q1+q2);
							if(q1>=q2){
								quality[i]=(byte) Tools.min(48, q1+1+q2/4);
							}else{
								quality[i]=(byte) Tools.min(48, q2+1+q1/4);
							}
						}
					}else if(b1=='N'){
						bases[i]=b2;
						quality[i]=q2;
					}else if(b2=='N'){
						//do nothing
					}else if(mergeVectors){
						if(q1<1 && q2<1){
							//Special case - e.g. Illumina calls bases at 0 quality.
							//Possibly best to keep the matching allele if one matches the ref.
							//But for now, do nothing.
							//This was causing problems changing perfect match strings into imperfect matches.
						}else if(q1==q2){
							assert(b1!=b2);
							bases[i]='N';
							quality[i]=0;
						}else if(q1>q2){
							bases[i]=b1;
							quality[i]=(byte)(q1-q2/2);
						}else{
							bases[i]=b2;
							quality[i]=(byte)(q2-q1/2);
						}
						assert(quality[i]>=0 && quality[i]<=48);
					}
				}
			}
		}

		//TODO:
		//Note that the read may need to be realigned after merging, so the match string may be rendered incorrect.

		if(mergeN && match!=null){
			if(r.match==null){match=null;}else{
				if(match.length!=r.match.length){match=null;}else{
					boolean ok=true;
					for(int i=0; i<match.length && ok; i++){
						byte a=match[i], b=r.match[i];
						if(a!=b){
							if((a=='m' || a=='S') && b=='N'){
								//do nothing;
							}else if(a=='N' && (b=='m' || b=='S')){
								match[i]=b;
							}else{
								ok=false;
							}
						}
					}
					if(!ok){match=null;}
				}
			}
		}

		if(mergeMate && mate!=null){
			mate.mergePrivate(r.mate, mergeVectors, mergeN, false);
			assert(copies==mate.copies);
		}
		assert(copies>1);

		assert(r!=this);
		assert(r!=this.mate);
		assert(r!=r.mate);
		assert(this!=this.mate);
		assert(r.mate==null || r.mate.mate==r);
		assert(this.mate==null || this.mate.mate==this);
		assert(r.mate==null || r.numericID==r.mate.numericID);
		assert(mate==null || numericID==mate.numericID);
	}

	/** Formats this read with the native-text serializer. */
	@Override
	public String toString(){return toText(false).toString();}

	/** Returns candidate sites as text in a new builder, or a dot when none are emitted. */
	public ByteBuilder toSites(){return toSites((ByteBuilder)null);}

	/** Appends candidate-site text without a final newline, allocating sb when null.
	 * Null sites are skipped but can leave extra tabs after the first emitted site.
	 * An empty/no-emitted-sites result appends a dot; existing builder contents remain. */
	public ByteBuilder toSites(ByteBuilder sb){
		if(numSites()==0){
			if(sb==null){sb=new ByteBuilder(2);}
			sb.append('.');
		}else{
			if(sb==null){sb=new ByteBuilder(sites.size()*20);}
			int appended=0;
			for(SiteScore ss : sites){
				if(appended>0){sb.append('\t');}
				if(ss!=null){
					ss.toBytes(sb);
					appended++;
				}
			}
			if(appended==0){sb.append('.');}
		}
		return sb;
	}

	/** Returns obj directly when its class is exactly ByteBuilder; otherwise builds its text or an empty builder.
	 * The direct-return case shares mutable storage with the caller's payload. */
	public ByteBuilder toInfo(){
		if(obj==null){return new ByteBuilder();}
		if(obj.getClass()==ByteBuilder.class){return (ByteBuilder)obj;}
		return new ByteBuilder(obj.toString());
	}

	/** Appends obj's text to bb, requiring a nonnull builder when obj is present.
	 * Unlike toInfo(), copies a ByteBuilder payload into bb; absent obj returns bb unchanged. */
	public ByteBuilder toInfo(ByteBuilder bb){
		if(obj==null){return bb;}
		if(obj.getClass()==ByteBuilder.class){return bb.append((ByteBuilder)obj);}
		return bb.append(obj.toString());
	}

	/** Formats through FASTQ.toFASTQ using current global output policy.
	 * @return New builder containing one record without its final quality-line newline */
	public ByteBuilder toFastq(){return FASTQ.toFASTQ(this, (ByteBuilder)null);}

	/** Appends a FASTQ record without a final newline; null bb creates a builder. */
	public ByteBuilder toFastq(ByteBuilder bb){return FASTQ.toFASTQ(this, bb);}

	/** Formats FASTA with Shared.FASTA_WRAP and no final newline.
	 * @return New builder containing the header and any sequence lines */
	public ByteBuilder toFasta(){return toFasta(Shared.FASTA_WRAP);}
	/** Appends FASTA with Shared.FASTA_WRAP; null bb creates a builder. */
	public ByteBuilder toFasta(ByteBuilder bb){return toFasta(Shared.FASTA_WRAP, bb);}

	/** Formats FASTA in a new builder, using the requested width and no final newline. */
	public ByteBuilder toFasta(int wrap){return toFasta(wrap, (ByteBuilder)null);}
	
	/** Appends id:bases:qualities with fixed ASCII-64 encoding and no newline.
	 * Requires bb. Missing qualities use FAKE_QUAL uniformly; does not use FASTQ output offset. */
	public ByteBuilder toScarf(ByteBuilder bb){
		bb.append(id).colon();
		bb.append(bases).colon();
		if(quality!=null){bb.appendQuality(quality, 64);}else{
			for(int i=0; i<length(); i++){bb.append((byte)(Shared.FAKE_QUAL+64));}
		}
		return bb;
	}

	/**
	 * Appends FASTA with the requested line width and no final newline.
	 * A null ID uses numericID; null/empty bases emit only the header. Existing builder
	 * contents remain. This does not serialize mates or mapping metadata.
	 * @param wrap Line wrap length (0 or negative for no wrapping)
	 * @param bb Optional ByteBuilder to append to (null creates new one)
	 * @return ByteBuilder containing FASTA representation
	 */
	public ByteBuilder toFasta(int wrap, ByteBuilder bb){
		if(wrap<1){wrap=Integer.MAX_VALUE;}
		int len=(id==null ? Tools.stringLength(numericID) : id.length())+(bases==null ? 0 : bases.length+bases.length/wrap)+5;
		if(bb==null){bb=new ByteBuilder(len+1);}
		bb.append('>');
		if(id==null){bb.append(numericID);}else{bb.append(id);}
		if(bases!=null){
			int pos=0;
			while(pos<bases.length-wrap){
				bb.append('\n');
				bb.append(bases, pos, wrap);
				pos+=wrap;
			}
			if(pos<bases.length){
				bb.append('\n');
				bb.append(bases, pos, bases.length-pos);
			}
		}
		return bb;
	}

	/** Builds a fresh SamLine from this Read and its pair number.
	 * @return New SAM text builder without a final newline */
	public ByteBuilder toSam(){return toSam((ByteBuilder)null);}

	/** Constructs a new SamLine rather than directly serializing the retained samline object.
	 * Appends to bb (allocating when null) without a final newline. */
	public ByteBuilder toSam(ByteBuilder bb){
		SamLine sl=new SamLine(this, pairnum());
		//		System.err.println("Called toSam on read "+id+"; num="+numericID+", pairnum="+pairnum()+"; result="+sl.toString());
		return sl.toBytes(bb);
	}

	/** Serializes native Read text in a new builder; see the mutating compression contract below. */
	public ByteBuilder toText(boolean okToCompressMatch){return toText(okToCompressMatch, (ByteBuilder)null);}

	/** Appends native text with fixed ASCII-33 qualities, gaps and optional sites; no final newline.
	 * Serializes selected fields, not obj, samline, mate links or all runtime metadata.
	 * When requested and enabled, temporarily replaces match and its representation bit,
	 * restoring them on normal return. Restoration is not protected by finally.
	 * A caller sharing this Read must exclude concurrent mutation/serialization. */
	public ByteBuilder toText(boolean okToCompressMatch, ByteBuilder bb){

		final byte[] oldmatch=match;
		final boolean oldshortmatch=this.shortmatch();
		if(COMPRESS_MATCH_BEFORE_WRITING && !shortmatch() && okToCompressMatch){
			match=toShortMatchString(match);
			setShortMatch(true);
		}

		if(bb==null){bb=new ByteBuilder();}
		bb.append(id);
		bb.tab();
		bb.append(numericID);
		bb.tab();
		bb.append(chrom);
		bb.tab();
		bb.append(Shared.strandCodes2[strand()]);
		bb.tab();
		bb.append(start);
		bb.tab();
		bb.append(stop);
		bb.tab();

		for(int i=maskArray.length-1; i>=0; i--){
			bb.append(flagToNumber(maskArray[i]));
		}
		bb.tab();

		bb.append(copies);
		bb.tab();

		bb.append(errors);
		bb.tab();
		bb.append(mapScore);
		bb.tab();

		if(bases==null){bb.append('.');}else{bb.append(bases);}
		bb.tab();

		//		int qualSum=0;
		//		int qualMin=99999;

		if(quality==null){
			bb.append('.');
		}else{
			bb.ensureExtra(quality.length);
			for(int i=0, j=bb.length; i<quality.length; i++, j++){
				byte q=quality[i];
				bb.array[j]=(byte)(q+ASCII_OFFSET);
				//				qualSum+=q;
				//				qualMin=Tools.min(q, qualMin);
			}
			bb.length+=quality.length;
		}
		bb.tab();

		if(insert<1){bb.append('.');}else{bb.append(insert);};
		bb.tab();

		if(true || quality==null){
			bb.append('.');
			bb.tab();
		}else{
			//			//These are not really necessary...
			//			sb.append(qualSum/quality.length);
			//			sb.append('\t');
		}

		if(match==null){bb.append('.');}else{bb.append(match);}
		bb.tab();

		if(gaps==null){
			bb.append('.');
		}else{
			for(int i=0; i<gaps.length; i++){
				if(i>0){bb.append('~');}
				bb.append(gaps[i]);
			}
		}

		if(sites!=null && sites.size()>0){

			assert(absdif(start, stop)<3000 || (gaps==null) == (sites.get(0).gaps==null)) :
				"\n"+this.numericID+"\n"+Arrays.toString(gaps)+"\n"+sites.toString()+"\n";

			for(SiteScore ss : sites){
				bb.tab();
				if(ss==null){
					bb.append((byte[])null);
				}else{
					ss.toBytes(bb);
				}
				//#008-fix [stream/Read#008]: removed a redundant second write here - `bb.append(ss==null?"null":ss.toText())` - that DUPLICATED the toBytes/null output just written above. SiteScore.toText and toBytes emit IDENTICAL content (verified), so every site was serialized TWICE into one column, breaking native-format round-trip (fromText's SiteScore.fromText would parse a doubled string). Kept the toBytes path, consistent with the originalSite block below. FLAG FOR BRIAN: validate the sited-read .bread round-trip; the ss==null element case is a pre-existing edge.
			}
		}

		if(originalSite!=null){
			bb.tab();
			bb.append('*');
			originalSite.toBytes(bb);
		}

		match=oldmatch;
		setShortMatch(oldshortmatch);

		return bb;
	}

	/** Extends uppercase-N runs of at least minGapIn to at least minGapOut, including terminal runs.
	 * Added bases receive Q0 when qualities exist. Replaces arrays only when length grows;
	 * leaves IDs, coordinates, match and mate metadata unchanged. Requires positive minGapIn. */
	public void inflateGaps(int minGapIn, int minGapOut){
		assert(minGapIn>0);
		if(!containsNocalls()){return;}
		final ByteBuilder bbb=new ByteBuilder();
		final ByteBuilder bbq=(quality==null ? null : new ByteBuilder());

		int gap=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			byte q=(quality==null ? 0 : quality[i]);
			if(b=='N'){
				gap++;
			}else{
				while(gap>=minGapIn && gap<minGapOut){
					gap++;
					bbb.append('N');
					if(bbq!=null){bbq.append(0);}
				}
				gap=0;
			}
			bbb.append(b);
			if(bbq!=null){bbq.append(q);}
		}

		while(gap>=minGapIn && gap<minGapOut){//Handle trailing bases
			gap++;
			bbb.append('N');
			if(bbq!=null){bbq.append(0);}
		}

		assert(bbb.length()>=bases.length);
		if(bbb.length()>bases.length){
			bases=bbb.toBytes();
			if(bbq!=null){quality=bbq.toBytes();}
		}
	}

	/** Splits nonnull bases on uppercase-N runs into fresh Read slices named id_c1, id_c2, etc.
	 * Each constructor receives copied bases/qualities and the original numericID, not mapping flags.
	 * Only slices meeting minContig enter the returned list, but numbering includes omitted slices.
	 * When agp is true, stores generated row bytes in this.obj, including rows for filtered contigs
	 * and internal/trailing gaps; leading gaps have no row. Otherwise leaves obj unchanged.
	 * Asserts that obj is initially empty; does not change this read's sequence arrays. */
	public ArrayList<Read> breakAtGaps(final boolean agp, final int minContig){
		ArrayList<Read> list=new ArrayList<Read>();
		byte prev='N';
		int lastN=-1, lastBase=-1;
		int contignum=1;
		long feature=1;
		ByteBuilder bb=(agp ? new ByteBuilder() : null);
		assert(obj==null);
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b=='N'){
				if(prev!='N'){
					final int start=lastN+1, stop=i;
					byte[] b2=KillSwitch.copyOfRange(bases, start, stop);
					byte[] q2=(quality==null ? null : KillSwitch.copyOfRange(quality, start, stop));
					Read r=new Read(b2, q2, id+"_c"+contignum, numericID);
					if(r.length()>=minContig){list.add(r);}
					contignum++;

					if(bb!=null){
						bb.append(id).append('\t');
						bb.append(start+1).append('\t');
						bb.append(stop).append('\t');
						bb.append(feature).append('\t');
						feature++;
						bb.append('W').append('\t');
						bb.append(r.id).append('\t');
						bb.append(1).append('\t');
						bb.append(r.length()).append('\t');
						bb.append('+').append('\n');
					}
				}
				lastN=i;
			}else{
				if(bb!=null && prev=='N' && lastBase>=0){
					bb.append(id).append('\t');
					bb.append(lastBase+2).append('\t');
					bb.append(i).append('\t');
					bb.append(feature).append('\t');
					feature++;
					bb.append('N').append('\t');
					bb.append((i-lastBase-1)).append('\t');
					bb.append("scaffold").append('\t');
					bb.append("yes").append('\t');
					bb.append("paired-ends").append('\n');
				}
				lastBase=i;
			}
			prev=b;
		}
		if(prev!='N'){
			final int start=lastN+1, stop=bases.length;
			byte[] b2=KillSwitch.copyOfRange(bases, start, stop);
			byte[] q2=(quality==null ? null : KillSwitch.copyOfRange(quality, start, stop));
			Read r=new Read(b2, q2, id+"_c"+contignum, numericID);
			if(r.length()>=minContig){list.add(r);}
			contignum++;

			if(bb!=null){
				bb.append(id).append('\t');
				bb.append(start+1).append('\t');
				bb.append(stop).append('\t');
				bb.append(feature).append('\t');
				feature++;
				bb.append('W').append('\t');
				bb.append(r.id).append('\t');
				bb.append(1).append('\t');
				bb.append(r.length()).append('\t');
				bb.append('+').append('\n');
			}
		}else{
			if(bb!=null && prev=='N' && lastBase>=0){
				bb.append(id).append('\t');
				bb.append(lastBase+2).append('\t');
				bb.append(bases.length).append('\t');
				bb.append(feature).append('\t');
				feature++;
				bb.append('N').append('\t');
				bb.append((bases.length-lastBase-1)).append('\t');
				bb.append("scaffold").append('\t');
				bb.append("yes").append('\t');
				bb.append("paired-ends").append('\n');
			}
			lastBase=bases.length;
		}
		if(bb!=null){obj=bb.toBytes();}
		return list;
	}

	/** Reverse-complements nonnull bases, reverses present qualities and toggles the strand bit.
	 * Uses the current Vector backend; leaves coordinates, match, mate and SamLine unchanged.
	 * @return This read after in-place mutation */
	public Read reverseComplement(){
		Vector.reverseComplementInPlace(bases);
		Vector.reverseInPlace(quality);
		setStrand(strand()^1);
		return this;
	}

	/** Uses Vector's fast reverse-complement path, which preserves supported IUPAC complements.
	 * Requires nonnull bases, reverses present qualities and toggles the strand bit; leaves
	 * other alignment/mate metadata unchanged. Currently shares the backend of reverseComplement().
	 * @return This read after in-place mutation */
	public Read reverseComplementFast(){
		Vector.reverseComplementInPlaceFast(bases);
		Vector.reverseInPlace(quality);
		setStrand(strand()^1);
		return this;
	}

	/** Complements bases in place with the extended alphabet table; null bases are a no-op.
	 * Does not reverse arrays or change qualities, strand, mapping fields or mate state. */
	public void complement(){AminoAcid.complementBasesInPlace(bases);}

	/** Compares chromosome, start, stop and strand by int subtraction, in that order.
	 * Other metadata and sequence contents are ignored; requires nonnull o. */
	@Override
	public int compareTo(Read o){
		if(chrom!=o.chrom){return chrom-o.chrom;}
		if(start!=o.start){return start-o.start;}
		if(stop!=o.stop){return stop-o.stop;}
		if(strand()!=o.strand()){return strand()-o.strand();}
		return 0;
	}

	/** Creates a SiteScore from current mapping state and assigns it to originalSite.
	 * Shares gaps and match arrays. Score initialization depends on the paired flag;
	 * this is not a side-effect-free snapshot or a deep copy. */
	public SiteScore toSite(){
		assert(start<=stop) : this.toText(false);
		SiteScore ss=new SiteScore(chrom, strand(), start, stop, 0, 0, rescued(), perfect());
		if(paired()){
			ss.setSlowPairedScore(mapScore-1, mapScore);
		}else{
			ss.setSlowPairedScore(mapScore, 0);
		}
		ss.setScore(mapScore);
		ss.gaps=gaps;
		ss.match=match;
		originalSite=ss;
		return ss;
	}

	/** Returns the first retained SiteScore, or null for a null/empty list; does not sort. */
	public SiteScore topSite(){
		final SiteScore ss=(sites==null || sites.isEmpty()) ? null : sites.get(0);
		assert(sites==null || sites.isEmpty() || ss!=null) : "Top site is null for read "+this;
		return ss;
	}

	/** Returns the current candidate-list size, or zero for a null list. */
	public int numSites(){return (sites==null ? 0 : sites.size());}

	/** Creates and retains a new original site through toSite, sharing its arrays. */
	public SiteScore makeOriginalSite(){
		originalSite=toSite();
		return originalSite;
	}

	/** Copies coordinates, strand, slow score, rescue/perfect flags and shared arrays from ss.
	 * Gap repair may replace both this.gaps and ss.gaps. Does not set mapped/paired bits,
	 * change originalSite, or synchronize a SamLine with the copied coordinates. */
	public void setFromSite(SiteScore ss){
		assert(ss!=null);
		chrom=ss.chrom;
		setStrand(ss.strand);
		start=ss.start;
		stop=ss.stop;
		mapScore=ss.slowScore;
		setRescued(ss.rescued);
		gaps=ss.gaps;
		setPerfect(ss.perfect);

		match=ss.match;

		if(gaps!=null){
			gaps=ss.gaps=GapTools.fixGaps(start, stop, gaps, Shared.MINGAP);
			//			gaps[0]=Tools.min(gaps[0], start);
			//			gaps[gaps.length-1]=Tools.max(gaps[gaps.length-1], stop);
		}
	}

	//	public static int[] fixGaps(int a, int b, int[] gaps, int minGap){
	////		System.err.println("fixGaps input: "+a+", "+b+", "+Arrays.toString(gaps)+", "+minGap);
	//		int[] r=GapTools.fixGaps(a, b, gaps, minGap);
	////		System.err.println("fixGaps output: "+Arrays.toString(r));
	//		return r;
	//	}

	/** Applies the retained original site; requires a nonnull originalSite. */
	public void setFromOriginalSite(){setFromSite(originalSite);}
	/** Applies the first site and marks mapped, or clears basic site state and mapped on no site.
	 * The no-site branch does not clear every field; see clearSite(). */
	public void setFromTopSite(){
		final SiteScore ss=topSite();
		if(ss==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);
		setFromSite(ss);
	}

	/** Selects among leading score/flag-compatible sites using numericID when ambiguity handling is enabled.
	 * Assumes the candidate list is ordered by score.
	 * May swap that site to index zero. The unfinished mate-guided path clears paired
	 * flags on both reads and retries the primary path; maxPairDist is currently unused. */
	public void setFromTopSite(boolean randomIfAmbiguous, boolean primary, int maxPairDist){
		final SiteScore ss0=topSite();
		if(ss0==null){
			clearSite();
			setMapped(false);
			return;
		}
		setMapped(true);

		if(sites.size()==1 || !randomIfAmbiguous || !ambiguous()){
			setFromSite(ss0);
			return;
		}

		if(primary || mate==null || !mate.mapped() || !mate.paired()){
			int count=1;
			for(int i=1; i<sites.size(); i++){
				SiteScore ss=sites.get(i);
				if(ss.score<ss0.score || (ss0.perfect && !ss.perfect) || (ss0.semiperfect && !ss.semiperfect)){break;}
				count++;
			}

			int x=(int)(numericID%count);
			if(x>0){
				SiteScore ss=sites.get(x);
				sites.set(0, ss);
				sites.set(x, ss0);
			}
			setFromSite(sites.get(0));
			return;
		}

		//		assert(false) : "TODO: Proper strand orientation, and more.";
		//TODO: Also, this code appears to sometimes duplicate sitescores(?)
		//		for(int i=0; i<list.size(); i++){
		//			SiteScore ss=list.get(i);
		//			if(ss.chrom==mate.chrom && Tools.min(Tools.absdifUnsigned(ss.start, mate.stop), Tools.absdifUnsigned(ss.stop, mate.start))<=maxPairDist){
		//				list.set(0, ss);
		//				list.set(i, ss0);
		//				setFromSite(ss);
		//				return;
		//			}
		//		}

		//If unsuccessful, recur unpaired.

		this.setPaired(false);
		mate.setPaired(false);
		setFromTopSite(randomIfAmbiguous, true, maxPairDist);
	}

	/** Clears mapping on this read and its linked mate, retaining the mate links. */
	public void clearPairMapping(){
		clearMapping();
		if(mate!=null){mate.clearMapping();}
	}

	/** Clears basic site state, match/candidate list and mapped/paired bits.
	 * Also clears the mate's paired bit, but retains mate links, originalSite and samline. */
	public void clearMapping(){
		clearSite();
		match=null;
		sites=null;
		setMapped(false);
		setPaired(false);
		if(mate!=null){mate.setPaired(false);}
	}

	/** Resets chrom/start/stop, strand, mapScore and gaps only.
	 * Does not clear mapped/paired/perfect bits, match, candidate sites, originalSite or samline. */
	public void clearSite(){
		chrom=-1;
		setStrand(0);
		start=-1;
		stop=-1;
		//		errors=0;
		mapScore=0;
		gaps=null;
	}


	/** Clears basic site state, match and candidates, retaining only synthetic/pair-number/swap bits.
	 * Optionally repeats this on the mate. Retains links, originalSite, samline and obj. */
	public void clearAnswers(boolean clearMate){
		//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
		clearSite();
		match=null;
		sites=null;
		flags=(flags&(SYNTHMASK|PAIRNUMMASK|SWAPMASK));
		if(clearMate && mate!=null){
			mate.clearSite();
			mate.match=null;
			mate.sites=null;
			mate.flags=(mate.flags&(SYNTHMASK|PAIRNUMMASK|SWAPMASK));
		}
		//		assert(mate==null || (pairnum()==0 && mate.pairnum()==1)) : pairnum()+", "+mate.pairnum();
	}


	/** Tests an unpaired alignment for incompatible mate coordinates/strands.
	 * Returns false without a mate, when paired() is already set, or if either read is unmapped.
	 * Different chromosomes or a start-minus-other-stop inner gap greater than maxdist are bad.
	 * requireCorrectStrands enforces the same/opposite choice in sameStrandPairs. Independently,
	 * opposite-strand mode rejects a plus-strand start at or beyond the minus-strand stop. */
	public boolean isBadPair(boolean requireCorrectStrands, boolean sameStrandPairs, int maxdist){
		if(mate==null || paired()){return false;}
		if(!mapped() || !mate.mapped()){return false;}
		if(chrom!=mate.chrom){return true;}

		{
			int inner;
			if(start<=mate.start){inner=mate.start-stop;}else{inner=start-mate.stop;}
			if(inner>maxdist){return true;}
		}
		//		if(absdif(start, mate.start)>maxdist){return true;}
		if(requireCorrectStrands){
			if((strand()==mate.strand())!=sameStrandPairs){return true;}
		}
		if(!sameStrandPairs){
			if(strand()==Shared.PLUS && mate.strand()==Shared.MINUS){
				if(start>=mate.stop){return true;}
			}else if(strand()==Shared.MINUS && mate.strand()==Shared.PLUS){
				if(mate.start>=stop){return true;}
			}
		}
		return false;
	}

	/** Counts literal uppercase S bytes in nonnull match; does not expand run counts. */
	public int countMismatches(){
		assert(match!=null);
		int x=0;
		for(byte b : match){
			if(b=='S'){x++;}
		}
		return x;
	}

	/** Sums nucleotide-defined windows across this read and its mate, if present.
	 * Each read uses numValidKmers's null, alphabet and nonpositive-k behavior.
	 * @param k Desired window length; normally positive
	 * @return Sum of the two per-read window counts */
	public int numValidPairKmers(int k){return numValidKmers(k)+(mate==null ? 0 : mate.numValidKmers(k));}

	/** Counts length-k windows containing only defined nucleotide symbols; ignores the amino flag.
	 * Null bases return zero; nonnull bytes must index the nucleotide lookup table.
	 * @param k Desired positive window length; nonpositive values count every array position
	 * @return Number of windows under this policy */
	public int numValidKmers(int k){
		if(bases==null){return 0;}
		int len=0, counted=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			if(x<0){len=0;}else{len++;}
			if(len>=k){counted++;}
		}
		return counted;
	}

	/** Calculates sequence identity based on alignment match string.
	 * @return Identity fraction from 0.0 to 1.0 */
	public final float identity(){return identity(match);}

	/** Tests this expanded match for a same-symbol I/X/Y run longer than maxlen. */
	public final boolean hasLongInsertion(int maxlen){return hasLongInsertion(match, maxlen);}

	/** Tests this expanded match for a D run longer than maxlen. */
	public final boolean hasLongDeletion(int maxlen){return hasLongDeletion(match, maxlen);}

	/** Accumulates encoded match lengths except D/C/X/Y/d; requires mapped state and arrays.
	 * Accepts expanded symbols or numeric run counts without changing the match array.
	 * Returns zero when unmapped, match/bases are null, or the match array is empty. */
	public int mappedNonClippedBases(){
		if(!mapped() || match==null || bases==null){return 0;}

		int len=0;
		//The initial sentinel is not a run; exclude it when the first symbol arrives.
		//Otherwise coverage callers receive an extra base for every nonempty match.
		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='0' || mode=='D' || mode=='C' || mode=='X' || mode=='Y' || mode=='d'){

					}else{
						len+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='D' || mode=='C' || mode=='X' || mode=='Y' || mode=='d'){

			}else{
				len+=current;
			}
			mode=c;
			current=0;
		}
		return len;
	}

	/** Tests the header's chastity indicator with malformed-field diagnostics enabled. */
	public boolean failsChastity(){return failsChastity(true);}

	/** Tests Y in the first space-delimited Illumina filter field, including the slash-prefixed form.
	 * Null, spaceless or too-short headers return false. Malformed available fields return false
	 * when processAssertions is false; otherwise assertions/fatal diagnostics apply. Does not change id. */
	public boolean failsChastity(boolean processAssertions){
		if(id==null){return false;}
		int space=id.indexOf(' ');
		if(space<0 || space+5>id.length()){return false;}
		char a=id.charAt(space+1);
		char b=id.charAt(space+2);
		char c=id.charAt(space+3);
		char d=id.charAt(space+4);

		if(a=='/'){
			if(b<'1' || b>'4' || c!=':'){
				if(!processAssertions){return false;}
				KillSwitch.kill("Strangely formatted read.  Please disable chastityfilter with the flag chastityfilter=f.  id:"+id);
			}
			return d=='Y';
		}else{
			if(processAssertions){
				assert(a=='1' || a=='2' || a=='3' || a=='4') : id;
				assert(b==':') : id;
				assert(d==':');
			}
			if(a<'1' || a>'4' || b!=':' || d!=':'){
				if(!processAssertions){return false;}
				KillSwitch.kill("Strangely formatted read.  Please disable chastityfilter with the flag chastityfilter=f.  id:"+id);
			}
			return c=='Y';
		}
	}

	/** Returns true when the whole suffix after the final colon fails a whitelist or base/plus syntax check.
	 * Null id returns false even when failIfNoBarcode is true. An absent/misplaced colon
	 * returns failIfNoBarcode; this method does not strip trailing comments like headerToBarcode. */
	public boolean failsBarcode(HashSet<String> set, boolean failIfNoBarcode){
		if(id==null){return false;}

		final int loc=id.lastIndexOf(':');
		final int loc2=Tools.max(id.indexOf(' '), id.indexOf('/'));
		if(loc<0 || loc<=loc2 || loc>=id.length()-1){
			return failIfNoBarcode;
		}

		if(set==null){
			for(int i=loc+1; i<id.length(); i++){
				char c=id.charAt(i);
				boolean ok=(c=='+' || AminoAcid.isFullyDefined(c));
				if(!ok){return true;}
			}
			return false;
		}else{
			String code=id.substring(loc+1);
			return !set.contains(code);
		}
	}

	/** Extracts this read's barcode, allocating a parser for nonnull id; see headerToBarcode for missing-header policy. */
	public String barcode(boolean failIfNoBarcode){return headerToBarcode(id, failIfNoBarcode, null);}

	/** Extracts this read's barcode, reusing/mutating nonnull ihp when id is present. */
	public String barcode(boolean failIfNoBarcode, IlluminaHeaderParser2 ihp){return headerToBarcode(id, failIfNoBarcode, ihp);}

	/**
	 * @return The rname of this Read's SamLine, if present and mapped.
	 */
	public String rnameS(){
		if(samline==null || !samline.mapped()){return null;}
		return samline.rnameS();
	}

	/** Dispatches to the configured probability or arithmetic-score average.
	 * countUndefined affects only the probability branch; maxBases below one means all qualities. */
	public double avgQuality(boolean countUndefined, int maxBases){return AVERAGE_QUALITY_BY_PROBABILITY ? avgQualityByProbabilityDouble(countUndefined, maxBases) : avgQualityByScoreDouble(maxBases);}

	/** Dispatches to the configured integer average: rounded probability-derived Phred
	 * versus truncating arithmetic mean. countUndefined affects only the probability branch. */
	public int avgQualityInt(boolean countUndefined, int maxBases){return AVERAGE_QUALITY_BY_PROBABILITY ? avgQualityByProbabilityInt(countUndefined, maxBases) : avgQualityByScoreInt(maxBases);}

	/** Returns zero for absent/empty bases; otherwise delegates to the static probability mean. */
	public int avgQualityByProbabilityInt(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return avgQualityByProbabilityInt(bases, quality, countUndefined, maxBases);
	}

	/** Returns zero for absent/empty bases; otherwise delegates to the static probability mean. */
	public double avgQualityByProbabilityDouble(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return avgQualityByProbabilityDouble(bases, quality, countUndefined, maxBases);
	}

	/** Multiplies modeled base-correctness probabilities; this is not an average quality.
	 * Returns zero for absent/empty bases before calling the static helper. */
	public double probabilityErrorFree(boolean countUndefined, int maxBases){
		if(bases==null || bases.length==0){return 0;}
		return probabilityErrorFree(bases, quality, countUndefined, maxBases);
	}

	/** Truncating arithmetic mean of selected qualities, flooring negative values to zero.
	 * Returns zero for absent/empty bases and 40 for null qualities. maxBases below one
	 * selects all qualities; a nonempty read requires at least one quality when present. */
	public int avgQualityByScoreInt(int maxBases){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null){return 40;}
		int x=0, limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		for(int i=0; i<limit; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/limit;
	}

	/** Arithmetic mean with double division, flooring negative qualities to zero.
	 * Returns zero for absent/empty bases and 40 for null qualities. maxBases below one
	 * selects all qualities. Present-but-empty qualities with nonempty bases yield NaN. */
	public double avgQualityByScoreDouble(int maxBases){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null){return 40;}
		int x=0, limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		for(int i=0; i<limit; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/(double)limit;
	}

	/** Truncating mean of the first n qualities, flooring negative values to zero.
	 * Empty bases or n exceeding quality length return zero; null qualities or n below
	 * one return 40 after the empty-base check. Used by BBMap tipsearch. */
	public int avgQualityFirstNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 40;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=0; i<n; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}

	/** Truncating mean of the last n qualities, using base length for the starting index.
	 * Requires matching base/quality lengths. Defaults match avgQualityFirstNBases. */
	public int avgQualityLastNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 40;}
		assert(quality!=null);
		int x=0;
		if(n>quality.length){return 0;}
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			x+=(b<0 ? 0 : b);
		}
		return x/n;
	}

	/** Minimum quality capped above at 41; absent bases/qualities or empty qualities return 41. */
	public int minQuality(){
		byte min=41;
		if(bases!=null && quality!=null){
			for(byte q : quality){
				min=Tools.min(min, q);
			}
		}
		return min;
	}

	/** Minimum among the first n qualities, without flooring negative values.
	 * Empty bases or n exceeding quality length return zero; null qualities or n below
	 * one return 41 after the empty-base check. Used by BBMap tipsearch. */
	public byte minQualityFirstNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 41;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[0];
		for(int i=1; i<n; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}

	/** Minimum among the last n qualities, indexed by base length; requires matching lengths.
	 * Defaults match minQualityFirstNBases. */
	public byte minQualityLastNBases(int n){
		if(bases==null || bases.length==0){return 0;}
		if(quality==null || n<1){return 41;}
		assert(quality!=null && n>0);
		if(n>quality.length){return 0;}
		byte x=quality[bases.length-n];
		for(int i=bases.length-n; i<bases.length; i++){
			byte b=quality[i];
			if(b<x){x=b;}
		}
		return x;
	}

	/** Tests for a byte above '9' other than m; asserts nonnull match and valid flag state. */
	public boolean containsNonM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m'){return true;}
		}
		return false;
	}

	/** Tests for a byte above '9' other than m/N, skipping encoded count digits. */
	public boolean containsNonNM(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N'){return true;}
		}
		return false;
	}

	/** Tests for a byte above '9' other than m/N/C; this is a symbol filter, not variant calling. */
	public boolean containsVariants(){
		assert(match!=null && valid()) : (match==null)+", "+(valid())+"\n"+samline+"\n";
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N' && b!='C'){return true;}
		}
		return false;
	}

	/** Tests for leading C or a final C run before trailing count digits; ignores internal clips. */
	public boolean containsClipping(){
		assert(match!=null && valid()) : (match==null)+", "+(valid())+"\n"+samline+"\n";
		if(match.length<1){return false;}
		if(match[0]=='C'){return true;}
		for(int i=match.length-1; i>0; i--){
			if(match[i]=='C'){return true;}
			if(match[i]>'9'){break;}
		}
		return false;
	}

	/** Counts m/S/I lengths from this match using the static encoded-token helper. */
	public int countAlignedBases(){return countAlignedBases(match);}

	//Ignores N/clip; counts each encoded D token once (D3 once, DDD three times).
	/**
	 * Counts alignment errors from match string.
	 * Counts substitution/insertion lengths and one per encoded deletion token.
	 * Compress consecutive D symbols first when a contiguous deletion should count once.
	 * @return Number of errors
	 */
	public int countErrors(){return countErrors(match);}

	/** Returns this match's m,S,C,N,I,D lengths in that order. */
	public int[] countMatchSymbols(){return countMatchSymbols(match);}

	/** Tests for a byte above '9' other than m/N/X/Y. */
	public boolean containsNonNMXY(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='N' && b!='X' && b!='Y'){return true;}
		}
		return false;
	}

	/** Tests for literal S/s/D/I symbols; does not count or decode run lengths. */
	public boolean containsSDI(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S' || b=='s' || b=='D' || b=='I'){return true;}
		}
		return false;
	}

	/** Tests for a byte above '9' other than m/s/N/S. */
	public boolean containsNonNMS(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b>'9' && b!='m' && b!='s' && b!='N' && b!='S'){return true;}
		}
		return false;
	}

	/** Tests an expanded match for at least num consecutive uppercase S symbols.
	 * Asserts nonnull/valid/non-short state; does not decode numeric runs. */
	public boolean containsConsecutiveS(int num){
		assert(match!=null && valid() && !shortmatch());
		int cnt=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			assert(b!='M');
			if(b=='S'){
				cnt++;
				if(cnt>=num){return true;}
			}else{
				cnt=0;
			}
		}
		return false;
	}

	/** Tests for literal I/D/X/Y symbols without interpreting their counts. */
	public boolean containsIndels(){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' || b=='D' || b=='X' || b=='Y'){return true;}
		}
		return false;
	}

	/** Counts uppercase S lengths in this expanded or run-counted match. */
	public int countSubs(){
		assert(match!=null && valid()) : (match!=null)+", "+valid()+", "+shortmatch();
		return countSubs(match);
	}

	/** Tests literal byte equality with c, including count digits when requested. */
	public boolean containsInMatch(char c){
		assert(match!=null && valid());
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b==c){return true;}
		}
		return false;
	}

	/** Tests nonnull bases for literal uppercase N, independently of the amino flag. */
	public boolean containsNocalls(){
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='N'){return true;}
		}
		return false;
	}

	/** Counts literal uppercase N in nonnull bases, independently of the amino flag. */
	public int countNocalls(){return countNocalls(bases);}

	/** Delegates to the nucleotide homopolymer helper; requires nonnull bases and ignores the amino flag. */
	public int longestHomopolymer(){return longestHomopolymer(bases);}

	/** Tests the nucleotide ACGTN lookup for rejected symbols; null bases return false.
	 * This ignores the amino flag and requires bytes that index the lookup table. */
	public boolean containsNonACGTN(){
		if(bases==null){return false;}
		for(byte b : bases){
			if(AminoAcid.baseToNumberACGTN[b]<0){return true;}
		}
		return false;
	}

	/** Tests for symbols undefined in this read's nucleotide or amino alphabet.
	 * Null bases return false; nonnull bytes must index the selected lookup table. */
	public boolean containsUndefined(){
		if(bases==null){return false;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		for(byte b : bases){
			if(symbolToNumber[b]<0){return true;}
		}
		return false;
	}

	/** Tests for ASCII a..z in bases; null bases return false. */
	public boolean containsLowercase(){
		if(bases==null){return false;}
		for(byte b : bases){
			if(Tools.isLowerCase(b)){return true;}
		}
		return false;
	}

	/** Counts symbols undefined in this read's nucleotide or amino alphabet.
	 * Null bases return zero; nonnull bytes must index the selected lookup table. */
	public int countUndefined(){
		if(bases==null){return 0;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		int n=0;
		for(byte b : bases){
			n+=(symbolToNumber[b]>=0 ? 0 : 1);
			//			if(symbolToNumber[b]<0){n++;}
		}
		return n;
	}

	/** Tests for at least min consecutive defined symbols in this read's alphabet.
	 * Null bases return min&lt;=0. For nonpositive min, nonnull arrays still need a defined symbol.
	 * Empty/all-undefined arrays return false; bytes must index the selected alphabet table. */
	public boolean hasMinConsecutiveBases(final int min){
		if(bases==null){return min<=0;}
		final byte[] symbolToNumber=AminoAcid.symbolToNumber(amino());
		int len=0;
		for(byte b : bases){
			if(symbolToNumber[b]<0){len=0;}else{
				len++;
				if(len>=min){return true;}
			}
		}
		return false;
	}


	/** Returns the smallest literal uppercase A/C/G/T count; other symbols are ignored.
	 * Null bases return zero. This does not change with the amino flag. */
	public int minBaseCount(){
		if(bases==null){return 0;}
		int a=0, c=0, g=0, t=0;
		for(byte b : bases){
			if(b=='A'){a++;}else if(b=='C'){c++;}else if(b=='G'){g++;}else if(b=='T'){t++;}
		}
		return Tools.min(a, c, g, t);
	}

	/** Tests this nonnull, valid-flagged match for a literal X or Y anywhere. */
	public boolean containsXY(){
		assert(match!=null && valid());
		return containsXY(match);
	}

	/** Checks only first-byte X or last-byte Y; asserts agreement with a full scan when valid. */
	public boolean containsXY2(){
		if(match==null || match.length<1){return false;}
		boolean b=(match[0]=='X' || match[match.length-1]=='Y');
		assert(!valid() || b==containsXY());
		return b;
	}

	/** Checks terminal X/Y as in containsXY2, or C at either literal endpoint.
	 * Does not skip trailing numeric counts when checking the last byte. */
	public boolean containsXYC(){
		if(match==null || match.length<1){return false;}
		boolean b=(match[0]=='X' || match[match.length-1]=='Y');
		assert(!valid() || b==containsXY());
		return b || match[0]=='C' || match[match.length-1]=='C';
	}

	/** Replaces 'B' in match string with 'S', 'm', or 'N' */
	public boolean fixMatchB(){
		assert(match!=null);
		final ChromosomeArray ca;
		if(Data.GENOME_BUILD>=0){
			ca=Data.getChromosome(chrom);
		}else{
			ca=null;
		}
		boolean originallyShort=shortmatch();
		if(originallyShort){match=toLongMatchString(match);}
		int mloc=0, cloc=0, rloc=start;
		for(; mloc<match.length; mloc++){
			byte m=match[mloc];

			if(m=='B'){
				byte r=(ca==null ? (byte)'?' : ca.get(rloc));
				byte c=bases[cloc];
				if(r=='N' || c=='N'){
					match[mloc]='N';
				}else if(r==c || Tools.toUpperCase(r)==Tools.toUpperCase(c)){
					match[mloc]='m';
				}else{
					if(ca==null){
						if(originallyShort){
							match=toShortMatchString(match);
						}
						for(int i=0; i<match.length; i++){
							if(match[i]=='B'){match[i]='N';}
						}
						return false;
					}
					match[mloc]='S';
				}
				cloc++;
				rloc++;
			}else if(m=='m' || m=='S' || m=='N' || m=='s' || m=='C'){
				cloc++;
				rloc++;
			}else if(m=='D'){
				rloc++;
			}else if(m=='I' || m=='X' || m=='Y'){
				cloc++;
			}
		}
		if(originallyShort){match=toShortMatchString(match);}
		return true;
	}

	/** Delegates to the reverse-end expected-error sum using this read's arrays. */
	public float expectedTipErrors(boolean countUndefined, int maxBases){return expectedTipErrors(bases, quality, countUndefined, maxBases);}

	/** Sums expected errors for this read and its linked mate; absent qualities contribute zero. */
	public float expectedErrorsIncludingMate(boolean countUndefined){
		float a=expectedErrors(countUndefined, length());
		float b=(mate==null ? 0 : mate.expectedErrors(countUndefined, mate.length()));
		return a+b;
	}

	/** Delegates to the forward expected-error sum using this read's arrays. */
	public float expectedErrors(boolean countUndefined, int maxBases){return expectedErrors(bases, quality, countUndefined, maxBases);}

	/** Counts X/Y symbols and low-quality S symbols in an expanded match, stopping at the match or base-array end.
	 * Returns zero without qualities; otherwise needs nonnull bases/match and qualities for visited bases.
	 * m/s/N/X/Y/I consume one base; D and unrecognized bytes do not. Each X/Y counts one error;
	 * S counts only below quality 19. No match expansion or probability weighting is performed. */
	public int estimateErrors(){
		if(quality==null){return 0;}
		assert(match!=null) : this.toText(false);

		int count=0;
		for(int ci=0, mi=0; ci<bases.length && mi<match.length; mi++){

			//			byte b=bases[ci];
			byte q=quality[ci];
			byte m=match[mi];
			if(m=='m' || m=='s' || m=='N'){
				ci++;
			}else if(m=='X' || m=='Y'){
				ci++;
				count++;
			}else if(m=='I'){
				ci++;
			}else if(m=='D'){

			}else if(m=='S'){
				ci++;
				if(q<19){
					count++;
				}
			}

		}
		return count;
	}

	/** Counts expanded m/S/D/I/N plus D runs reaching max(minSplice,1).
	 * C joins N and X/Y join I. If numeric counts are encountered with shortmatch set,
	 * expands match in place, clears the flag, prints a warning and retries.
	 * @return New counts in m,S,D,I,N-or-C,qualifying-D-runs order */
	public int[] countErrors(int minSplice){
		assert(match!=null) : this.toText(false);
		int m=0;
		int s=0;
		int d=0;
		int i=0;
		int n=0;
		int splice=0;

		byte prev=' ';
		int streak=0;
		minSplice=Tools.max(minSplice, 1);

		for(int pos=0; pos<match.length; pos++){
			final byte b=match[pos];

			if(b==prev){streak++;}else{streak=1;}

			if(b=='m'){
				m++;
			}else if(b=='N' || b=='C'){
				n++;
			}else if(b=='X' || b=='Y'){
				i++;
			}else if(b=='I'){
				i++;
			}else if(b=='D'){
				d++;
				if(streak==minSplice){splice++;}
			}else if(b=='S'){
				s++;
			}else{
				if(Tools.isDigit(b) && shortmatch()){
					System.err.println("Warning! Found read in shortmatch form during countErrors():\n"+this); //Usually caused by verbose output.
					if(mate!=null){System.err.println("mate:\n"+mate.id+"\t"+new String(mate.bases));}
					System.err.println("Stack trace: ");
					new Exception().printStackTrace();
					match=toLongMatchString(match);
					setShortMatch(false);
					return countErrors(minSplice);
				}else{
					throw new RuntimeException("\nUnknown symbol "+(char)b+":\n"+new String(match)+"\n"+this+"\nshortmatch="+this.shortmatch());
				}
			}

			prev=b;
		}

		//		assert(i==0) : i+"\n"+this+"\n"+new String(match)+"\n"+Arrays.toString(new int[] {m, s, d, i, n, splice});

		return new int[] {m, s, d, i, n, splice};
	}

	/** Compresses this match and sets SHORTMATCHMASK unless already marked short.
	 * An already-short input asserts when doAssertion is true, otherwise returns unchanged. */
	public void toShortMatchString(boolean doAssertion){
		if(shortmatch()){
			assert(!doAssertion);
			return;
		}
		match=toShortMatchString(match);
		setShortMatch(true);
	}

	/** Expands this match and clears SHORTMATCHMASK unless already marked long.
	 * An already-long input asserts when doAssertion is true, otherwise returns unchanged. */
	public void toLongMatchString(boolean doAssertion){
		if(!shortmatch()){
			assert(!doAssertion);
			return;
		}
		match=toLongMatchString(match);
		setShortMatch(false);
	}

	/** Returns the reference name parsed from this nonnull SYN identifier using pairnum().
	 * CustomHeader's caught parse exceptions may disable FASTQ.PARSE_CUSTOM and leave the name null;
	 * this helper asserts the SYN prefix but is not a complete header validator. */
	public String parseCustomRname(){
		assert(id.startsWith("SYN")) : "Can't parse header "+id;
		return new CustomHeader(id, pairnum()).rname;
	}

	/** Returns obj cast to FloatList, without copying or checking its runtime type first. */
	public FloatList fetchVector(){
		FloatList x=(FloatList)obj;
		return x;
	}

	/** Returns obj cast to ByteBuilder, sharing its mutable contents. */
	public ByteBuilder fetchBB(){
		ByteBuilder x=(ByteBuilder)obj;
		return x;
	}

	/** Returns the caller-owned auxiliary object without transferring or copying it. */
	public Object obj(){return obj;}

	/** Assigns and returns x; asserts that x is nonnull and the old slot is empty. */
	public Object setObj(Object x){
		assert(x!=null);
		Object old=obj;
		assert(old==null) : obj.getClass()+", "+x.getClass(); //Just a warning.  This should probably not happen.
		obj=x;
		return obj;
	}

	/** Clears the auxiliary slot and returns its previous value. */
	public Object nullifyObj(){
		Object old=obj;
		obj=null;
		return old;
	}

	/** Clears the auxiliary slot and returns its previous value. */
	public Object nullifyObject(){
		Object old=obj;
		obj=null;
		return old;
	}

	/** Alias for nullifyObject(); returns the previous auxiliary value. */
	public Object setObjNull(){return nullifyObject();}

	/** Returns obj cast to TrimRead; no copy or type conversion is performed. */
	public TrimRead fetchTrimRead(){return (TrimRead)obj;}

	/** Returns obj's runtime class, or null for an empty slot. */
	public Class getObjectClass(){return obj==null ? null : obj.getClass();}

	/** Returns the same auxiliary reference as obj(). */
	public Object fetchObject(){return obj;}

	/** Returns the stored identifier, which may be null. */
	public String name(){return id;}

	/** Returns obj as a Long timestamp; asserts that the payload is a nonnull Long. */
	public long time(){
		assert(obj!=null && obj.getClass()==Long.class) : obj;
		return ((Long)obj).longValue();
	}
	/** Number of bases in this pair, including the mate if present. */
	public int pairLength(){return length()+mateLength();}
	/** Number of length-k windows across this read and its mate, without base-validity filtering. */
	public int numPairKmers(int k){return numKmers(k)+numMateKmers(k);}
	/** Number of reads in this pair.  Returns 1 if the read has no mate, and 2 if it does. */
	public int pairCount(){return 1+mateCount();}
	/** Counts stored mapped bits on this read and its linked mate, from zero to two. */
	public int pairMappedCount(){return (mapped() ? 1 : 0)+(mate==null || !mate.mapped() ? 0 : 1);}
	/** Returns sequence-array length, or zero when bases are null. */
	public int length(){return bases==null ? 0 : bases.length;}
	/** Counts length-k windows without base-validity filtering; null bases yield zero. */
	public int numKmers(int k){return bases==null ? 0 : Tools.max(0, bases.length-k+1);}
	/** Returns quality-array length, or zero when qualities are absent. */
	public int qlength(){return quality==null ? 0 : quality.length;}
	/** Returns the mate's sequence length, or zero without a mate. */
	public int mateLength(){return mate==null ? 0 : mate.length();}
	/** Counts the mate's length-k windows without validity filtering, or zero without a mate. */
	public int numMateKmers(int k){return mate==null ? 0 : mate.numKmers(k);}
	/** Returns the mate's identifier, or null without a mate or mate identifier. */
	public String mateId(){return mate==null ? null : mate.id;}
	/** Returns one if the mate reference is nonnull, otherwise zero. */
	public int mateCount(){return mate==null ? 0 : 1;}
	/** Tests the linked mate's stored mapped bit; false without a mate. */
	public boolean mateMapped(){return mate==null ? false : mate.mapped();}
	/** Tests the stored mapped bit on either this read or its linked mate. */
	public boolean eitherMapped(){return mapped() || mateMapped();}
	/** Returns the mate's approximate memory size, or zero without a mate. */
	public long countMateBytes(){return mate==null ? 0 : mate.countBytes();}
	/** Returns the mate's field-based FASTQ size estimate, or zero without a mate. */
	public long countMateFastqBytes(){return mate==null ? 0 : mate.countFastqBytes();}
	/** Number of bytes this read pair uses in memory, approximately */
	public long countPairBytes(){return countBytes()+(mate==null ? 0 : mate.countBytes());}

	/** Estimates memory using fixed object overhead plus selected arrays, ID and SAM data.
	 * Does not traverse mate, sites or gaps, or account for shared references exactly. */
	public long countBytes(){
		long sum=144; //Approximate per-read overhead
		sum+=(bases==null ? 0 : bases.length+16);
		sum+=(quality==null ? 0 : quality.length+16);
		sum+=(id==null ? 0 : id.length()*2+16);
		sum+=(match==null ? 0 : match.length+16);
		sum+=(samline==null ? 0 : samline.countBytes());
		sum+=(obj==null ? 0 : 32);
		return sum;
	}

	/** Estimates FASTQ bytes as six framing bytes plus existing bases, qualities and ID lengths.
	 * Does not account for writer options, generated qualities/IDs or character encoding. */
	public long countFastqBytes(){
		long sum=6;//4 newlines, +, @
		sum+=(bases==null ? 0 : bases.length);
		sum+=(quality==null ? 0 : quality.length);
		sum+=(id==null ? 0 : id.length());
		return sum;
	}

	/** Counts exact leading bytes after narrowing base to byte; requires nonnull bases. */
	public int countLeading(final char base){return countLeft((byte)base);}
	/** Counts exact trailing bytes after narrowing base to byte; requires nonnull bases. */
	public int countTrailing(final char base){return countRight((byte)base);}
	/** Alias for countLeft(byte) after narrowing base to byte. */
	public int countLeft(final char base){return countLeft((byte)base);}
	/** Alias for countRight(byte) after narrowing base to byte. */
	public int countRight(final char base){return countRight((byte)base);}

	/** Counts the leading run equal to base, without case folding; requires nonnull bases. */
	public int countLeft(final byte base){
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b!=base){return i;}
		}
		return bases.length;
	}

	/** Counts the trailing run equal to base, without case folding; requires nonnull bases. */
	public int countRight(final byte base){
		for(int i=bases.length-1; i>=0; i--){
			final byte b=bases[i];
			if(b!=base){return bases.length-i-1;}
		}
		return bases.length;
	}

	/** Dispatches an exact TrimRead payload's untrim operation and clears obj.
	 * Returns true for that payload type even if the operation reports no restoration;
	 * returns false and leaves obj unchanged for any other payload. */
	public boolean untrim(){
		if(obj==null || obj.getClass()!=TrimRead.class){return false;}
		((TrimRead)obj).untrim();
		obj=null;
		return true;
	}

	/** Counts trailing ASCII a-z bytes; requires nonnull bases. */
	public int trailingLowerCase(){
		for(int i=bases.length-1; i>=0;){
			if(Tools.isLowerCase(bases[i])){
				i--;
			}else{
				return bases.length-i-1;
			}
		}
		return bases.length;
	}
	/** Counts leading ASCII a-z bytes; requires nonnull bases. */
	public int leadingLowerCase(){
		for(int i=0; i<bases.length; i++){
			if(!Tools.isLowerCase(bases[i])){return i;}
		}
		return bases.length;
	}

	/** Looks up the stored strand bit in Shared.strandCodes2. */
	public char strandChar(){return Shared.strandCodes2[strand()];}
	/** Returns strand: 0 for plus strand, 1 for minus strand */
	public byte strand(){return (byte)(flags&1);}
	/** Tests the stored mapped bit; does not validate coordinates or alignment data. */
	public boolean mapped(){return (flags&MAPPEDMASK)==MAPPEDMASK;}
	/** Returns the stored paired-alignment bit, independently of whether mate is nonnull. */
	public boolean paired(){return (flags&PAIREDMASK)==PAIREDMASK;}
	/** Tests the stored synthetic-read bit. */
	public boolean synthetic(){return (flags&SYNTHMASK)==SYNTHMASK;}
	/** Tests the stored ambiguous-mapping bit. */
	public boolean ambiguous(){return (flags&AMBIMASK)==AMBIMASK;}
	/** Tests the stored perfect-alignment bit without examining match or bases. */
	public boolean perfect(){return (flags&PERFECTMASK)==PERFECTMASK;}
	//	public boolean semiperfect(){return perfect() ? true : list!=null && list.size()>0 ? list.get(0).semiperfect : false;} //TODO: This is a hack.  Add a semiperfect flag.
	/** Tests the stored rescued-alignment bit. */
	public boolean rescued(){return (flags&RESCUEDMASK)==RESCUEDMASK;}
	/** Tests the stored discard bit without examining the mate. */
	public boolean discarded(){return (flags&DISCARDMASK)==DISCARDMASK;}
	/** Tests the stored invalid bit without performing validation. */
	public boolean invalid(){return (flags&INVALIDMASK)==INVALIDMASK;}
	/** Tests the stored swapped bit without comparing arrays or mate state. */
	public boolean swapped(){return (flags&SWAPMASK)==SWAPMASK;}
	/** Tests the short-match representation bit without inspecting match bytes. */
	public boolean shortmatch(){return (flags&SHORTMATCHMASK)==SHORTMATCHMASK;}
	/** Tests the insert-valid bit without checking the stored insert value. */
	public boolean insertvalid(){return (flags&INSERTMASK)==INSERTMASK;}
	/** Tests the stored adapter bit. */
	public boolean hasAdapter(){return (flags&ADAPTERMASK)==ADAPTERMASK;}
	/** Tests the stored secondary-alignment bit. */
	public boolean secondary(){return (flags&SECONDARYMASK)==SECONDARYMASK;}
	/** Tests the stored supplementary-alignment bit. */
	public boolean supplementary(){return (flags&SUPPLEMENTARYMASK)==SUPPLEMENTARYMASK;}
	/** True iff neither secondary (SAM 0x100) nor supplementary (SAM 0x800) is set; i.e.
	 * this is the SAM representative (primary) line of its read. */
	public boolean primary(){return (flags&(SECONDARYMASK|SUPPLEMENTARYMASK))==0;}
	/** Tests the amino-acid representation bit; does not inspect sequence symbols. */
	public boolean aminoacid(){return (flags&AAMASK)==AAMASK;}
	/** Alias for aminoacid(). */
	public boolean amino(){return (flags&AAMASK)==AAMASK;}
	/** Tests the stored junk bit, independently of the invalid/discard bits. */
	public boolean junk(){return (flags&JUNKMASK)==JUNKMASK;}
	/** Tests the validation-completed bit; junk/discard/invalid are independent flags. */
	public boolean validated(){return (flags&VALIDATEDMASK)==VALIDATEDMASK;}
	/** Tests the caller-managed tested bit. */
	public boolean tested(){return (flags&TESTEDMASK)==TESTEDMASK;}
	/** Tests the stored inverted-repeat bit. */
	public boolean invertedRepeat(){return (flags&IRMASK)==IRMASK;}
	/** Tests the stored trimmed bit without inspecting obj or sequence length. */
	public boolean trimmed(){return (flags&TRIMMEDMASK)==TRIMMEDMASK;}

	/** For paired ends: 0 for read1, 1 for read2 */
	public int pairnum(){return (flags&PAIRNUMMASK)>>PAIRNUMSHIFT;}
	/** Tests only that INVALIDMASK is clear; does not run validation or inspect junk. */
	public boolean valid(){return !invalid();}

	/** Tests whether all bits in mask are set; an empty mask returns true. */
	public boolean getFlag(int mask){return (flags&mask)==mask;}
	/** Returns one when all mask bits are set, otherwise zero. */
	public int flagToNumber(int mask){return (flags&mask)==mask ? 1 : 0;}

	/** Sets or clears all mask bits, preserving other flags and associated fields. */
	public void setFlag(int mask, boolean b){
		flags=(flags&~mask);
		if(b){flags|=mask;}
	}

	/** Sets the strand bit to zero or one (asserted), preserving other flags. */
	public void setStrand(int b){
		assert(b==1 || b==0);
		flags=(flags&(~1))|b;
	}

	/** For paired ends: 0 for read1, 1 for read2 */
	public void setPairnum(int b){
		//		System.err.println("Setting pairnum to "+b+" for "+id);
		//		assert(!id.equals("2_chr1_0_1853883_1853982_1845883_ecoli_K12") || b==1);
		assert(b==1 || b==0);
		flags=(flags&(~PAIRNUMMASK))|(b<<PAIRNUMSHIFT);
		//		assert(pairnum()==b);
	}

	/** Changes only the paired-alignment bit; does not link/unlink mates or modify the mate. */
	public void setPaired(boolean b){
		flags=(flags&~PAIREDMASK);
		if(b){flags|=PAIREDMASK;}
	}

	/** Changes only the synthetic bit; does not populate originalSite. */
	public void setSynthetic(boolean b){
		flags=(flags&~SYNTHMASK);
		if(b){flags|=SYNTHMASK;}
	}

	/** Changes only the ambiguous bit; does not modify candidate sites. */
	public void setAmbiguous(boolean b){
		flags=(flags&~AMBIMASK);
		if(b){flags|=AMBIMASK;}
	}

	/** Updates and returns the perfect bit using the top site's score/flag and match data.
	 * Clears the bit without a top site. A top score equal to maxScore or a perfect site
	 * sets it, with assertion-only consistency checks; otherwise tests the match directly. */
	public boolean setPerfectFlag(int maxScore){
		final SiteScore ss=topSite();
		if(ss==null){
			setPerfect(false);
		}else{
			assert(ss.slowScore<=maxScore) : maxScore+", "+ss.slowScore+", "+ss.toText();

			if(ss.slowScore==maxScore || ss.perfect){
				assert(testMatchPerfection(true)) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
				setPerfect(true);
			}else{
				boolean flag=testMatchPerfection(false);
				setPerfect(flag);
				assert(flag || !ss.perfect) : "flag="+flag+", ss.perfect="+ss.perfect+"\nmatch="+new String(match)+"\n"+this.toText(false);
				assert(!flag || ss.slowScore>=maxScore) : "\n"+ss+"\n"+maxScore+"\n"+this+"\n"+mate+"\n";
			}
		}
		return perfect();
	}

	/** Tests existing match representation and absence of uppercase N in nonnull bases.
	 * Null match returns the supplied default. Expanded form requires one m per base;
	 * short form is empty or starts with m and contains only m/digits; decoded length is unchecked. */
	private boolean testMatchPerfection(boolean returnIfNoMatch){
		if(match==null){return returnIfNoMatch;}
		boolean flag=(match.length==bases.length);
		if(shortmatch()){
			flag=(match.length==0 || match[0]=='m');
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m' || Tools.isDigit(match[i]));}
		}else{
			for(int i=0; i<match.length && flag; i++){flag=(match[i]=='m');}
		}
		for(int i=0; i<bases.length && flag; i++){flag=(bases[i]!='N');}
		return flag;
	}

	/** Returns GC/(AT+GC) over symbols recognized by baseToNumber, ignoring other symbols.
	 * Returns zero for null/empty bases or no recognized G/C; this is a nucleotide statistic. */
	public float gc(){
		if(bases==null || bases.length<1){return 0;}
		int at=0, gc=0;
		for(byte b : bases){
			int x=AminoAcid.baseToNumber[b];
			if(x>-1){
				if(x==0 || x==3){at++;}else{gc++;}
			}
		}
		if(gc<1){return 0;}
		return gc*1f/(at+gc);
	}

	/** Replaces exact swapFrom bytes in place and returns the number of matching positions.
	 * Returns zero for null bases; equal from/to bytes still count. Leaves qualities alone. */
	public int swapBase(byte swapFrom, byte swapTo){
		if(bases==null){return 0;}
		int swaps=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]==swapFrom){
				bases[i]=swapTo;
				swaps++;
			}
		}
		return swaps;
	}

	/** Replaces each base with remap[base] in place; null bases are a no-op.
	 * Requires every original byte to be a valid table index; leaves qualities alone. */
	public void remap(byte[] remap){
		if(bases==null){return;}
		for(int i=0; i<bases.length; i++){
			bases[i]=remap[bases[i]];
		}
	}

	/** Remaps bases in place and counts positions whose byte actually changes.
	 * Null bases return zero; each original byte must index remap. Qualities are unchanged. */
	public int remapAndCount(byte[] remap){
		if(bases==null){return 0;}
		int swaps=0;
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=remap[a];
			if(a!=b){
				bases[i]=b;
				swaps++;
			}
		}
		return swaps;
	}

	/** Replaces positions rejected by baseToNumberACGTN and zeros their qualities if present.
	 * A negative replacement selects every position without consulting that table.
	 * Otherwise each base must be a valid table index; qualities must cover selected positions.
	 * Returns selected-position count (even for equal replacements), or zero for null bases. */
	public int convertUndefinedTo(byte b){
		if(bases==null){return 0;}
		int changed=0;
		for(int i=0; i<bases.length; i++){
			if(b<0 || AminoAcid.baseToNumberACGTN[bases[i]]<0){
				changed++;
				bases[i]=b;
				if(quality!=null){quality[i]=0;}
			}
		}
		return changed;
	}

	/** Exchanges base and quality array references with the mate, leaving all metadata alone.
	 * Missing mate asserts and otherwise returns without changing the read. */
	public void swapBasesWithMate(){
		if(mate==null){
			assert(false);
			return;
		}
		byte[] temp=bases;
		bases=mate.bases;
		mate.bases=temp;
		temp=quality;
		quality=mate.quality;
		mate.quality=temp;
	}

	/** Returns the stored insert when its valid bit is set, otherwise -1. */
	public int insert(){return insertvalid() ? insert : -1;}

	/** Estimates insert size from this read and its linked mate using the mapped dispatcher. */
	public int insertSizeMapped(boolean ignoreStrand){return insertSizeMapped(this, mate, ignoreStrand);}

	/** Uses originalSite coordinates rather than the current mapping.
	 * Without a mate returns this site's inclusive span, or -1 without a site.
	 * With a mate returns zero if either site is absent, otherwise calls insertSize;
	 * this read must then have nonnull bases. Asserts symmetry when this pair number is below the mate's. */
	public int insertSizeOriginalSite(){
		if(mate==null){
			//			System.err.println("A: "+(originalSite==null ? "null" : (originalSite.stop-originalSite.start+1)));
			return (originalSite==null ? -1 : originalSite.stop-originalSite.start+1);
		}

		final SiteScore ssa=originalSite, ssb=mate.originalSite;
		final int x;
		if(ssa==null || ssb==null){
			//			System.err.println("B: 0");
			x=0;
		}else{
			x=insertSize(ssa, ssb, bases.length, mate.length());
		}

		assert(pairnum()>=mate.pairnum() || x==mate.insertSizeOriginalSite());
		return x;
	}

	/** Clones this read with copied base/quality slices [from,to) and no mate link.
	 * Other arrays/objects and mapping metadata retain shallow-copy semantics; coordinates
	 * and pairing flags are not adjusted for the slice. */
	public Read subRead(int from, int to){
		Read r=this.clone();
		r.bases=KillSwitch.copyOfRange(bases, from, to);
		r.quality=(quality==null ? null : KillSwitch.copyOfRange(quality, from, to));
		r.mate=null;
		//		assert(Tools.indexOf(r.bases, (byte)'J')<0);
		return r;
	}

	/** Joins this read and its already-oriented mate using the valid stored positive insert.
	 * Returns this unchanged when the insert/mate is unavailable. Asserts the historical
	 * small-insert plausibility check before delegating to the static joining method. */
	public Read joinRead(){
		if(insert<1 || mate==null || !insertvalid()){return this;}
		assert(insert>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, insert);
	}

	/** Joins this read and its already-oriented mate using x, independent of the insert-valid bit.
	 * Returns this unchanged for a missing mate or nonpositive x. */
	public Read joinRead(int x){
		if(x<1 || mate==null){return this;}
		//assert(x>9 || bases.length<20) : "Perhaps old read format is being used?  This appears to be a quality value, not an insert.\n"+this+"\n\n"+mate+"\n";
		return joinRead(this, mate, x);
	}

	/** Returns null below minlen, or a one-element list containing this read when no split is needed.
	 * The multipart branch is unfinished and asserts under -ea; its known quality-slice concern
	 * is annotated below. This method is not a validated general read-splitting operation.
	 * @param minlen Minimum accepted sequence length
	 * @param maxlen Positive target maximum fragment length
	 * @return Null, a list retaining this read, or the legacy multipart result when assertions are disabled */
	public ArrayList<Read> split(int minlen, int maxlen){
		int len=bases==null ? 0 : bases.length;
		if(len<minlen){return null;}
		int parts=(len+maxlen-1)/maxlen;
		ArrayList<Read> subreads=new ArrayList<Read>(parts);
		if(len<=maxlen){
			subreads.add(this);
		}else{
			float ideal=Tools.max(minlen, len/(float)parts);
			int actual=(int)ideal;
			//TODO: Possible bug [stream/Read#007] - the multi-part split path (len>maxlen, the PRIMARY purpose of split) is UNFINISHED: this active assert(false) crashes under -ea (always on) before any part is produced. Also ~L3914 subquals copyOfRange(quality, a, b+1) is off-by-one vs subbases (a,b) (mismatched length + possible AIOOBE). Both in this fenced path. FLAG FOR BRIAN (is split's multi-part path needed? trace callers).
			assert(false) : "TODO"; //Some assertion goes here, I forget what
			for(int i=0; i<parts; i++){
				int a=i*actual;
				int b=(i+1)*actual;
				if(b>bases.length){b=bases.length;}
				//				if(b-a<)
				byte[] subbases=KillSwitch.copyOfRange(bases, a, b);
				byte[] subquals=(quality==null ? null : KillSwitch.copyOfRange(quality, a, b+1));
				Read r=new Read(subbases, subquals, id+"_"+i, numericID, flags);
				subreads.add(r);
			}
		}
		return subreads;
	}

	/** Produces one value per length-k window, leaving -1 where a window contains an undefined base.
	 * For positive k at most 31, emits forward two-bit encodings or the maximum of forward/reverse
	 * encodings when makeCanonical is true. Larger k delegates to toLongKmers and its scratch rules.
	 * Requires gap=0 (positive values throw); null/short bases return null. Reuses kmers only when its length
	 * exactly matches the window count, overwriting every element. Does not change this read. */
	public long[] toKmers(final int k, final int gap, long[] kmers, boolean makeCanonical, Kmer longkmer){
		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
		if(k>31){return toLongKmers(k, kmers, makeCanonical, longkmer);}
		if(bases==null || bases.length<k+gap){return null;}

		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);

		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;

		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				kmers[i-k+1]=makeCanonical ? Tools.max(kmer, rkmer) : kmer;
			}
		}
		return kmers;
	}

	//	/** Generate and return an array of canonical kmers for this read */
	//	public long[] toKmers(final int k, final int gap, long[] kmers, boolean makeCanonical, Kmer longkmer) {
	//		if(gap>0){throw new RuntimeException("Gapped reads: TODO");}
	//		if(k>31){return toLongKmers(k, kmers, makeCanonical, longkmer);}
	//		if(bases==null || bases.length<k+gap){return null;}
	//
	//		final int kbits=2*k;
	//		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
	//
	//		int len=0;
	//		long kmer=0;
	//		final int arraylen=bases.length-k+1;
	//		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
	//		Arrays.fill(kmers, -1);
	//
	//		for(int i=0; i<bases.length; i++){
	//			byte b=bases[i];
	//			int x=AminoAcid.baseToNumber[b];
	//			if(x<0){
	//				len=0;
	//				kmer=0;
	//			}else{
	//				kmer=((kmer<<2)|x)&mask;
	//				len++;
	//
	//				if(len>=k){
	//					kmers[i-k+1]=kmer;
	//				}
	//			}
	//		}
	//
	////		System.out.println(new String(bases));
	////		System.out.println(Arrays.toString(kmers));
	//
	//		if(makeCanonical){
	//			this.reverseComplement();
	//			len=0;
	//			kmer=0;
	//			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
	//				byte b=bases[i];
	//				int x=AminoAcid.baseToNumber[b];
	//				if(x<0){
	//					len=0;
	//					kmer=0;
	//				}else{
	//					kmer=((kmer<<2)|x)&mask;
	//					len++;
	//
	//					if(len>=k){
	//						assert(kmer==AminoAcid.reverseComplementBinaryFast(kmers[j], k));
	//						kmers[j]=Tools.max(kmers[j], kmer);
	//					}
	//				}
	//			}
	//			this.reverseComplement();
	//
	////			System.out.println(Arrays.toString(kmers));
	//		}
	//
	//
	//		return kmers;
	//	}

	/** Produces Kmer.xor hashes per length-k window using caller-owned mutable scratch.
	 * Requires k greater than 31 and makeCanonical=true (asserted), with scratch configured for k.
	 * Null/short bases return null before touching scratch. Otherwise clears and advances scratch,
	 * reuses a correctly sized output array, and leaves undefined-base windows at -1.
	 * Values are hashes under Kmer's configuration, not reversible packed sequence encodings. */
	public long[] toLongKmers(final int k, long[] kmers, boolean makeCanonical, Kmer kmer){
		assert(k>31) : k;
		assert(makeCanonical);
		if(bases==null || bases.length<k){return null;}
		kmer.clear();

		final int arraylen=bases.length-k+1;
		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
		Arrays.fill(kmers, -1);

		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			kmer.addRight(b);
			if(!AminoAcid.isFullyDefined(b)){kmer.clear();}
			if(kmer.len>=k){
				kmers[i-k+1]=kmer.xor();
			}
		}

		return kmers;
	}

	/** Changes only the perfect bit; does not check or alter the alignment. */
	public void setPerfect(boolean b){
		flags=(flags&~PERFECTMASK);
		if(b){flags|=PERFECTMASK;}
	}

	/** Changes only the rescued bit. */
	public void setRescued(boolean b){
		flags=(flags&~RESCUEDMASK);
		if(b){flags|=RESCUEDMASK;}
	}

	/** Changes only the mapped bit; does not populate or clear alignment fields. */
	public void setMapped(boolean b){
		flags=(flags&~MAPPEDMASK);
		if(b){flags|=MAPPEDMASK;}
	}

	/** Changes the discard bit on this read and its linked mate, when present. */
	public void setPairDiscarded(boolean b){
		setDiscarded(b);
		if(mate!=null){mate.setDiscarded(b);}
	}
	/** Changes only this read's discard bit; leaves the mate untouched. */
	public void setDiscarded(boolean b){
		flags=(flags&~DISCARDMASK);
		if(b){flags|=DISCARDMASK;}
	}

	/** Changes only the invalid bit; does not run validation. */
	public void setInvalid(boolean b){
		flags=(flags&~INVALIDMASK);
		if(b){flags|=INVALIDMASK;}
	}

	/** Changes only the swapped bit; does not exchange arrays. */
	public void setSwapped(boolean b){
		flags=(flags&~SWAPMASK);
		if(b){flags|=SWAPMASK;}
	}

	/** Changes only the representation bit; does not encode or decode match. */
	public void setShortMatch(boolean b){
		flags=(flags&~SHORTMATCHMASK);
		if(b){flags|=SHORTMATCHMASK;}
	}

	/** Changes only the insert-valid bit; leaves the stored insert unchanged. */
	public void setInsertValid(boolean b){
		flags=(flags&~INSERTMASK);
		if(b){flags|=INSERTMASK;}
	}

	/** Changes only the adapter bit. */
	public void setHasAdapter(boolean b){
		flags=(flags&~ADAPTERMASK);
		if(b){flags|=ADAPTERMASK;}
	}

	/** Changes only the secondary bit; does not update the retained SamLine. */
	public void setSecondary(boolean b){
		flags=(flags&~SECONDARYMASK);
		if(b){flags|=SECONDARYMASK;}
	}

	/** Changes only the supplementary bit; does not update the retained SamLine. */
	public void setSupplementary(boolean b){
		flags=(flags&~SUPPLEMENTARYMASK);
		if(b){flags|=SUPPLEMENTARYMASK;}
	}

	/** Changes only the amino-acid bit; does not translate bases or alter qualities. */
	public void setAminoAcid(boolean b){
		flags=(flags&~AAMASK);
		if(b){flags|=AAMASK;}
	}

	/** Changes only the junk bit; does not alter invalid/discard state. */
	public void setJunk(boolean b){
		flags=(flags&~JUNKMASK);
		if(b){flags|=JUNKMASK;}
	}

	/** Changes only the validation-completed bit; does not run validation. */
	public void setValidated(boolean b){
		flags=(flags&~VALIDATEDMASK);
		if(b){flags|=VALIDATEDMASK;}
	}

	/** Changes only the caller-managed tested bit. */
	public void setTested(boolean b){
		flags=(flags&~TESTEDMASK);
		if(b){flags|=TESTEDMASK;}
	}

	/** Changes only the inverted-repeat bit. */
	public void setInvertedRepeat(boolean b){
		flags=(flags&~IRMASK);
		if(b){flags|=IRMASK;}
	}

	/** Changes only the trimmed bit; does not alter arrays or the auxiliary payload. */
	public void setTrimmed(boolean b){
		flags=(flags&~TRIMMEDMASK);
		if(b){flags|=TRIMMEDMASK;}
	}

	/** Stores a positive insert and marks it valid, or stores -1 and clears validity.
	 * Applies the same value/bit to the linked mate when present; leaves coordinates alone. */
	public void setInsert(int x){
		if(x<1){x=-1;}
		//		assert(x==-1 || x>9 || length()<20) : x+", "+length(); //Invalid assertion for synthetic reads.
		insert=x;
		setInsertValid(x>0);
		if(mate!=null){
			mate.insert=x;
			mate.setInsertValid(x>0);
		}
	}

	/** Returns the borrowed scaffold name selected by this mapped read's midpoint.
	 * Returns null when unmapped or when the optional single-scaffold check fails.
	 * Requires compatible loaded Data scaffold metadata; uses its padding-aware span policy.
	 * Does not copy the name or guarantee that both endpoints lie inside scaffold bases. */
	public byte[] getScaffoldName(boolean requireSingleScaffold){
		byte[] name=null;
		if(mapped()){
			if(!requireSingleScaffold || Data.isSingleScaffold(chrom, start, stop)){
				//TODO: Possible midpoint overflow - start+stop can overflow before division;
				//supported coordinate reachability needs checking, as in ScaffoldCoordinates.
				int idx=Data.scaffoldIndex(chrom, (start+stop)/2);
				name=Data.scaffoldNames[chrom][idx];
				//				int scaflen=Data.scaffoldLengths[chrom][idx];
				//				a1=Data.scaffoldRelativeLoc(chrom, start, idx);
				//				b1=a1-start1+stop1;
			}
		}
		return name;
	}

	/** Applies selected substitutions in place, classifying each original byte once.
	 * Uses baseToNumber (including lowercase and U as T), emits uppercase replacements,
	 * and leaves unselected/undefined symbols unchanged. Requires nonnull, table-indexable
	 * bases; qualities and alignment metadata are unchanged. */
	public void bisulfite(boolean AtoG, boolean CtoT, boolean GtoA, boolean TtoC){
		for(int i=0; i<bases.length; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x==0 && AtoG){bases[i]='G';}else if(x==1 && CtoT){bases[i]='T';}else if(x==2 && GtoA){bases[i]='A';}else if(x==3 && TtoC){bases[i]='C';}
		}
	}

	/** Copies bases/qualities/match/gaps and detaches the copied read from its mate.
	 * Clones the candidate list and SiteScore objects, but SiteScore.clone shares its arrays.
	 * originalSite is likewise shallow-cloned; obj and samline remain shared. Flags,
	 * including pairing state, are copied unchanged. This is not a complete deep copy. */
	public Read copy(){
		Read r=clone();
		r.bases=(r.bases==null ? null : r.bases.clone());
		r.quality=(r.quality==null ? null : r.quality.clone());
		r.match=(r.match==null ? null : r.match.clone());
		r.gaps=(r.gaps==null ? null : r.gaps.clone());
		r.originalSite=(r.originalSite==null ? null : r.originalSite.clone());
		r.sites=(ArrayList<SiteScore>) (r.sites==null ? null : r.sites.clone());
		r.mate=null;

		if(r.sites!=null){
			for(int i=0; i<r.sites.size(); i++){
				r.sites.set(i, r.sites.get(i).clone());
			}
		}
		return r;
	}

	/** Shallow clone: all arrays, mate/site links and auxiliary objects remain shared. */
	@Override
	public Read clone(){
		try{
			return (Read) super.clone();
		}catch(CloneNotSupportedException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}

	/** Shallow-clones an amino-flagged read and replaces its bases with AminoAcid.toNTs output.
	 * Clears the clone's amino bit and expands each quality to three copies of min(q+5,max quality).
	 * Other references, including mate and alignment state, retain clone sharing; coordinates
	 * and mate links are not rewritten for nucleotide space. The original read is unchanged.
	 * @return New read with converted bases and, when present, a new expanded quality array */
	public Read aminoToNucleic(){
		assert(aminoacid()) : "This read is not flagged as an amino acid sequence.";
		Read r=this.clone();
		r.setAminoAcid(false);
		r.bases=AminoAcid.toNTs(r.bases);
		if(quality!=null){
			byte[] ntquals=new byte[r.quality.length*3];
			for(int i=0; i<quality.length; i++){
				byte q=quality[i];
				byte q2=(byte)Tools.min(q+5, MAX_CALLED_QUALITY);
				final int j=i*3; ntquals[j]=ntquals[j+1]=ntquals[j+2]=q2;//#002-fix [stream/Read#002]: was ntquals[i]/[i+1]/[i+2] (index i not i*3) -> overlapping writes left 2/3 of ntquals zero + corrupted the rest; amino i maps to nt 3i,3i+1,3i+2
			}
			r.quality=ntquals;
		}
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns larger minus smaller using int arithmetic; does not prevent subtraction overflow. */
	private static final int absdif(int a, int b){return a>b ? a-b : b-a;}

	/** Returns historical native-text labels, not an exact schema of current toText output.
	 * The retained labels include a length column that the current serializer does not emit. */
	public static CharSequence header(){

		StringBuilder sb=new StringBuilder();
		sb.append("id");
		sb.append('\t');
		sb.append("numericID");
		sb.append('\t');
		sb.append("chrom");
		sb.append('\t');
		sb.append("strand");
		sb.append('\t');
		sb.append("start");
		sb.append('\t');
		sb.append("stop");
		sb.append('\t');

		sb.append("flags");
		sb.append('\t');

		sb.append("copies");
		sb.append('\t');

		sb.append("errors,fixed");
		sb.append('\t');
		sb.append("mapScore");
		sb.append('\t');
		sb.append("length");
		sb.append('\t');

		sb.append("bases");
		sb.append('\t');
		sb.append("quality");
		sb.append('\t');

		sb.append("insert");
		sb.append('\t');
		{
			//These are not really necessary...
			sb.append("avgQual");
			sb.append('\t');
		}

		sb.append("match");
		sb.append('\t');
		sb.append("SiteScores: "+SiteScore.header());
		return sb;
	}

	/**
	 * Parses the historical native-text layout using fixed field indexes.
	 * Builds new sequence/quality arrays, converts qualities from fixed ASCII-33 and
	 * invokes configured constructor validation. The strand comes from flag bits;
	 * selected legacy fields are ignored. Does not restore mate links, obj or samline.
	 * Site-column selection has the separate compatibility concern below. Optionally
	 * expands a match flagged as short. Null/empty input is not the null sentinel.
	 * @param line Nonnull tab-delimited native record, or the literal dot sentinel
	 * @return Parsed Read object, or null exactly when line equals "."
	 */
	public static Read fromText(String line){
		if(line.length()==1 && line.charAt(0)=='.'){return null;}

		String[] split=line.split("\t");

		//		if(split.length<17){
		//			throw new RuntimeException("Error parsing read from text.\n\n" +
		//					"This may be caused be attempting to parse the wrong format.\n" +
		//					"Please ensure that the file extension is correct:\n" +
		//					"\tFASTQ should end in .fastq or .fq\n" +
		//					"\tFASTA should end in .fasta or .fa, .fas, .fna, .ffn, .frn, .seq, .fsa\n" +
		//					"\tSAM should end in .sam\n" +
		//					"\tNative format should end in .txt or .bread\n" +
		//					"If a file is compressed, there must be a compression extension after the format extension:\n" +
		//					"\tgzipped files should end in .gz or .gzip\n" +
		//					"\tzipped files should end in .zip and have only 1 file per archive\n" +
		//					"\tbz2 files should end in .bz2\n");
		//		}

		final String id=new String(split[0]);
		long numericID=Long.parseLong(split[1]);
		//#004-fix [stream/Read#004]: was Byte.parseByte - chrom is an int and toText writes the full int (bb.append(chrom)); a native-format read with chrom>127 (multi-scaffold reference) threw NumberFormatException on round-trip. Now Integer.parseInt to match toText. Native (.bread/.txt) format; latent.
		int chrom=Integer.parseInt(split[2]);
		//		byte strand=Byte.parseByte(split[3]);
		int start=Integer.parseInt(split[4]);
		int stop=Integer.parseInt(split[5]);

		int flags=Integer.parseInt(split[6], 2);

		int copies=Integer.parseInt(split[7]);

		int errors;
		int errorsCorrected;
		if(split[8].indexOf(',')>=0){
			String[] estring=split[8].split(",");
			errors=Integer.parseInt(estring[0]);
			errorsCorrected=Integer.parseInt(estring[1]);
		}else{
			errors=Integer.parseInt(split[8]);
			errorsCorrected=0;
		}

		int mapScore=Integer.parseInt(split[9]);

		byte[] basesOriginal=split[10].getBytes();
		byte[] qualityOriginal=(split[11].equals(".") ? null : split[11].getBytes());

		if(qualityOriginal!=null){
			for(int i=0; i<qualityOriginal.length; i++){
				byte b=qualityOriginal[i];
				b=(byte) (b-ASCII_OFFSET);
				assert(b>=-1) : b;
				qualityOriginal[i]=b;
			}
		}

		int insert=-1;
		if(!split[12].equals(".")){insert=Integer.parseInt(split[12]);}

		byte[] match=null;
		if(!split[14].equals(".")){match=split[14].getBytes();}
		int[] gaps=null;
		if(!split[15].equals(".")){
			//#003-fix [stream/Read#003]: was split[16] - the gaps column is 15 (checked just above + written at col 15 by toText); split[16] is the first SiteScore, so a native-format gapped read parsed gaps from the wrong column (NumberFormatException / garbage). Now split[15].
			String[] gstring=split[15].split("~");
			gaps=new int[gstring.length];
			for(int i=0; i<gstring.length; i++){
				gaps[i]=Integer.parseInt(gstring[i]);
			}
		}

		//		assert(false) : split[16];

		Read r=new Read(basesOriginal, qualityOriginal, id, numericID, flags, chrom, start, stop);
		r.match=match;
		r.errors=errors;
		r.mapScore=mapScore;
		r.copies=copies;
		r.gaps=gaps;
		r.insert=insert;

		//TODO: Probable native-layout mismatch - toText starts appended site data at
		//column 16, but this decoder starts at 17/18, potentially skipping a site or
		//originalSite. Verify supported historical layouts before changing the selector.
		int firstScore=(ADD_BEST_SITE_TO_LIST_FROM_TEXT) ? 17 : 18;

		int scores=split.length-firstScore;

		int mSites=0;
		for(int i=firstScore; i<split.length; i++){
			if(split[i].charAt(0)!='*'){mSites++;}
		}

		//This can be disabled to handle very old text format.
		if(mSites>0){r.sites=new ArrayList<SiteScore>(mSites);}
		for(int i=firstScore; i<split.length; i++){
			SiteScore ss=SiteScore.fromText(split[i]);
			if(split[i].charAt(0)=='*'){r.originalSite=ss;}else{r.sites.add(ss);}
		}

		if(DECOMPRESS_MATCH_ON_LOAD && r.shortmatch()){
			r.toLongMatchString(true);
		}

		assert(r.numSites()==0 || absdif(r.start, r.stop)<3000 || (r.gaps==null) == (r.topSite().gaps==null)) :
			"\n"+r.numericID+", "+r.chrom+", "+r.strand()+", "+r.start+", "+r.stop+", "+Arrays.toString(r.gaps)+"\n"+r.sites+"\n"+line+"\n";

		return r;
	}

	/** Counts expanded symbols or numeric runs, grouping C/X/Y together and N/R together.
	 * Unrecognized symbol modes contribute nothing; this is not alphabet validation.
	 * @param match Expanded or run-counted match bytes
	 * @return Counts in m,S,D,I,C-or-X-or-Y,N-or-R order; null for null/empty input */
	public static final int[] matchToMsdicn(byte[] match){
		if(match==null || match.length<1){return null;}
		int[] msdicn=KillSwitch.allocInt1D(6);

		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						msdicn[0]+=current;
					}else if(mode=='S'){
						msdicn[1]+=current;
					}else if(mode=='D'){
						msdicn[2]+=current;
					}else if(mode=='I'){
						msdicn[3]+=current;
					}else if(mode=='C' || mode=='X' || mode=='Y'){
						msdicn[4]+=current;
					}else if(mode=='N' || mode=='R'){
						msdicn[5]+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				msdicn[0]+=current;
			}else if(mode=='S'){
				msdicn[1]+=current;
			}else if(mode=='D'){
				msdicn[2]+=current;
			}else if(mode=='I'){
				msdicn[3]+=current;
			}else if(mode=='C' || mode=='X' || mode=='Y'){
				msdicn[4]+=current;
			}else if(mode=='N' || mode=='R'){
				msdicn[5]+=current;
			}
		}
		return msdicn;
	}


	/** Calculates the reference-span convention used by alignment helpers.
	 * Counts m/S/D/C/X/N/R, ignoring I/Y and unrecognized modes. Supports expanded
	 * symbols and numeric runs. A final X run triggers the existing assertion.
	 * @param match Expanded or run-counted match bytes
	 * @return Calculated span, or zero for null/empty input */
	public static final int calcMatchLength(byte[] match){
		if(match==null || match.length<1){return 0;}

		byte mode='0', c='0';
		int current=0;
		int len=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						len+=current;
					}else if(mode=='S'){
						len+=current;
					}else if(mode=='D'){
						len+=current;
					}else if(mode=='I'){ //Do nothing
						//len+=current;
					}else if(mode=='C'){
						len+=current;
					}else if(mode=='X'){ //Not sure about this case, but adding seems fine
						len+=current;
						//						assert(false) : new String(match);
					}else if(mode=='Y'){ //Do nothing
						//len+=current;
						//						assert(false) : new String(match);
					}else if(mode=='N' || mode=='R'){
						len+=current;
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				len+=current;
			}else if(mode=='S'){
				len+=current;
			}else if(mode=='D'){
				len+=current;
			}else if(mode=='I'){ //Do nothing
				//len+=current;
			}else if(mode=='C'){
				len+=current;
			}else if(mode=='X'){ //Not sure about this case, but adding seems fine
				len+=current;
				assert(false) : new String(match);
			}else if(mode=='Y'){ //Do nothing
				//len+=current;
				//				assert(false) : new String(match);
			}else if(mode=='N' || mode=='R'){
				len+=current;
			}
		}
		return len;
	}

	/**
	 * Calculates sequence identity from match string.
	 * Uses FLAT_IDENTITY to select flat scoring or square-root deletion weighting;
	 * both branches count N/R with fractional weights. Does not modify match.
	 * @param match Expanded or run-counted alignment symbols
	 * @return Identity fraction from 0.0 to 1.0
	 */
	public static final float identity(byte[] match){
		if(FLAT_IDENTITY){
			return identityFlat(match, true);
		}else{
			return identitySkewed(match, true, true, false, false);
		}
	}

	/** Tests expanded bytes for an I, X or Y run strictly longer than maxlen.
	 * Different insertion-like symbols restart the run; numeric counts are not expanded.
	 * Null or an array shorter than maxlen returns false. */
	public static final boolean hasLongInsertion(byte[] match, int maxlen){
		if(match==null || match.length<maxlen){return false;}
		byte prev='0';
		int len=0;
		for(byte b : match){
			if(b=='I' || b=='X' || b=='Y'){
				if(b==prev){len++;}else{len=1;}
				if(len>maxlen){return true;}
			}else{
				len=0;
			}
			prev=b;
		}
		return false;
	}

	/** Tests expanded bytes for a D run strictly longer than maxlen; does not decode counts.
	 * Null or an array shorter than maxlen returns false. */
	public static final boolean hasLongDeletion(byte[] match, int maxlen){
		if(match==null || match.length<maxlen){return false;}
		byte prev='0';
		int len=0;
		for(byte b : match){
			if(b=='D'){
				if(b==prev){len++;}else{len=1;}
				if(len>maxlen){return true;}
			}else{
				len=0;
			}
			prev=b;
		}
		return false;
	}

	/** Computes flat identity from expanded or run-counted symbols.
	 * m contributes good length; S/I/X/Y/i/d and short D runs contribute bad length.
	 * D runs at least SamLine.INTRON_LIMIT and C/V are omitted. N/R contribute
	 * 0.25 good and 0.75 bad per base when penalized, otherwise nothing.
	 * Final-mode assertions are narrower than internal-mode assertions; this is not
	 * a general validator for arbitrary terminal symbols.
	 * @param match Match bytes, or null for zero
	 * @param penalizeN Include fractional N/R contributions
	 * @return Good weight divided by max(total weight,1); zero for null/empty input */
	public static final float identityFlat(byte[] match, boolean penalizeN){
		//		assert(false) : new String(match);
		if(match==null || match.length<1){return 0;}

		int good=0, bad=0, n=0;

		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						good+=current;
						//						System.out.println("G: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}else if(mode=='R' || mode=='N'){
						n+=current;
					}else if(mode=='C' || mode=='V'){
						//Do nothing
						//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
					}else if(mode!='0'){
						assert(mode=='S' || mode=='D' || mode=='I' || mode=='X' || mode=='Y' || mode=='i' || mode=='d') : (char)mode;
						if(mode!='D' || current<SamLine.INTRON_LIMIT){
							bad+=current;
						}
						//						System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				good+=current;
			}else if(mode=='R' || mode=='N'){
				n+=current;
			}else if(mode=='C' || mode=='V'){
				//Do nothing
				//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
			}else if(mode!='0'){
				assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
				if(mode!='D' || current<SamLine.INTRON_LIMIT){
					bad+=current;
				}
				//				System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
			}
		}


		float good2=good+n*(penalizeN ? 0.25f : 0);
		float bad2=bad+n*(penalizeN ? 0.75f : 0);
		float r=good2/Tools.max(good2+bad2, 1);
		//		assert(false) : new String(match)+"\nmode='"+(char)mode+"', current="+current+", good="+good+", bad="+bad;

		//		System.out.println("match="+new String(match)+"\ngood="+good+", bad="+bad+", r="+r);
		//		System.out.println(Arrays.toString(matchToMsdicn(match)));

		return r;
	}

	/** Computes identity with reduced internal D-run penalties and flat other errors.
	 * Expanded and run-counted forms are supported. Internal D runs below INTRON_LIMIT
	 * use ceil(sqrt(length)), ceil(log2(length)) or one, capped at their length;
	 * longer D runs and C/V are omitted. N/R use the same weights as identityFlat.
	 * The final-mode branch retains flat accumulation and narrower assertions.
	 * @param match Match bytes, or null for zero
	 * @param penalizeN Include 0.25 good/0.75 bad per N/R base
	 * @param sqrt Use square-root weighting before considering log
	 * @param log Use log2 weighting when sqrt is false
	 * @param single Unused; weight one is the fallback when both sqrt and log are false
	 * @return Good weight divided by max(total weight,1); zero for null/empty input */
	public static final float identitySkewed(byte[] match, boolean penalizeN, boolean sqrt, boolean log, boolean single){
		//		assert(false) : new String(match);
		if(match==null || match.length<1){return 0;}

		int good=0, bad=0, n=0;

		byte mode='0', c='0';
		int current=0;
		for(int i=0; i<match.length; i++){
			c=match[i];
			if(Tools.isDigit(c)){
				current=(current*10)+(c-'0');
			}else{
				if(mode==c){
					current=Tools.max(current+1, 2);
				}else{
					current=Tools.max(current, 1);

					if(mode=='m'){
						good+=current;
						//						System.out.println("G: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}else if(mode=='D'){
						if(current<SamLine.INTRON_LIMIT){
							int x;

							if(sqrt){x=(int)Math.ceil(Math.sqrt(current));}else if(log){x=(int)Math.ceil(Tools.log2(current));}else{x=1;}

							bad+=(Tools.min(x, current));
						}

						//						System.out.println("D: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad+", x="+x);
					}else if(mode=='R' || mode=='N'){
						n+=current;
					}else if(mode=='C' || mode=='V'){
						//Do nothing
						//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
					}else if(mode!='0'){
						assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y' || mode=='i' || mode=='d') : (char)mode;
						bad+=current;
						//						System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
					}
					mode=c;
					current=0;
				}
			}
		}
		if(current>0 || !Tools.isDigit(c)){
			current=Tools.max(current, 1);
			if(mode=='m'){
				good+=current;
			}else if(mode=='R' || mode=='N'){
				n+=current;
			}else if(mode=='C' || mode=='V'){
				//Do nothing
				//I assume this is clipped because it went off the end of a scaffold, and thus is irrelevant to identity
			}else if(mode!='0'){
				assert(mode=='S' || mode=='I' || mode=='X' || mode=='Y') : (char)mode;
				if(mode!='D' || current<SamLine.INTRON_LIMIT){
					bad+=current;
				}
				//				System.out.println("B: mode="+(char)mode+", c="+(char)c+", current="+current+", good="+good+", bad="+bad);
			}
		}


		float good2=good+n*(penalizeN ? 0.25f : 0);
		float bad2=bad+n*(penalizeN ? 0.75f : 0);
		float r=good2/Tools.max(good2+bad2, 1);
		//		assert(false) : new String(match)+"\nmode='"+(char)mode+"', current="+current+", good="+good+", bad="+bad;

		//		System.out.println("match="+new String(match)+"\ngood="+good+", bad="+bad+", r="+r);
		//		System.out.println(Arrays.toString(matchToMsdicn(match)));

		return r;
	}

	/** Uses IlluminaHeaderParser2 first, then a final-colon fallback for older header shapes.
	 * The fallback requires the colon after the first space/slash markers and cuts a trailing
	 * comment at the first following space, or at a tab only if there is no following space.
	 * Does not validate barcode symbols. Missing ID/barcode returns null unless the requested
	 * fatal policy is active with Shared.EA(). A supplied parser is overwritten only for nonnull id.
	 * @param id Header to parse; may be null
	 * @param failIfNoBarcode Request a fatal missing-barcode diagnostic when assertions are enabled
	 * @param ihp Reusable parser, or null to allocate one for a nonnull header
	 * @return Parsed barcode or null when absent under the nonfatal policy */
	public static String headerToBarcode(String id, boolean failIfNoBarcode, IlluminaHeaderParser2 ihp){

		if(id==null){
			if(failIfNoBarcode && Shared.EA()){
				KillSwitch.kill("Encountered a read header without a barcode:\n"+id+"\n");
			}
			return null;
		}

		if(ihp==null){ihp=new IlluminaHeaderParser2();}
		ihp.parse(id);
		String barcode=ihp.barcode();
		if(barcode!=null){return barcode;}
		
		final int loc=id.lastIndexOf(':');
		final int loc2=Tools.max(id.indexOf(' '), id.indexOf('/'));
		if(loc<0 || loc<=loc2 || loc>=id.length()-1){
			if(failIfNoBarcode && Shared.EA()){
				KillSwitch.kill("Encountered a read header without a barcode:\n"+id+"\n");
			}
			return null;
		}

		//This section allows for comments after the barcode
		final int bcStart=loc+1;
		int bcStop=id.indexOf(' ', bcStart);
		bcStop=(bcStop>=0 ? bcStop : id.indexOf('\t', bcStart));
		String code;
		if(bcStop<0){
			code=id.substring(bcStart);
		}else{
			code=id.substring(bcStart, bcStop);
		}

		return code;
	}

	/** Converts mean modeled error probability to rounded, capped integer Phred.
	 * Null qualities return 40; empty qualities return zero. The denominator is the
	 * selected quality count, even when undefined bases are excluded from the error sum.
	 * maxBases below one selects all qualities; arrays must contain matching positions. */
	public static int avgQualityByProbabilityInt(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 40;}
		if(quality.length==0){return 0;}
		float e=expectedErrors(bases, quality, countUndefined, maxBases);
		final int div=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		float p=e/div;
		return QualityTools.probErrorToPhred(p);
	}

	/** Converts mean modeled error probability to Phred through QualityTools' double conversion.
	 * The sum/mean intermediates are floats. Null qualities return 40, empty qualities
	 * return zero; the denominator includes skipped undefined positions. maxBases below
	 * one selects all qualities. The conversion saturates at 60 for very small probabilities. */
	public static double avgQualityByProbabilityDouble(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 40;}
		if(quality.length==0){return 0;}
		float e=expectedErrors(bases, quality, countUndefined, maxBases);
		final int div=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		float p=e/div;
		return QualityTools.probErrorToPhredDouble(p);
	}

	/** Sums m/S/I lengths in expanded or run-counted form; C/N/D contribute zero.
	 * Requires nonnull input and asserts the supported alphabet for completed tokens. */
	public static int countAlignedBases(byte[] match){
		//Short version; correct but slow
		//		int[] sym=countMatchSymbols(match);
		//		return sym[0]+sym[1]+sym[4];

		//Long version; faster
		int msi=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m' || last=='S' || last=='I'){
					msi+=current;
				}else{
					//Ignore
					assert(last=='?' || last=='C' || last=='N' || last=='D') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m' || last=='S' || last=='I'){
			msi+=current;
		}else{
			//Ignore
			assert(last=='?' || last=='C' || last=='N' || last=='D') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return msi;
	}

	/** Counts S/I lengths plus one per encoded D token, ignoring m/N/C.
	 * Requires nonnull input; compress repeated D first for contiguous deletion-event counts. */
	public static int countErrors(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m+=current;
				}else if(last=='S'){
					S+=current;
				}else if(last=='C'){
					C+=current;
				}else if(last=='N'){
					N+=current;
				}else if(last=='I'){
					I+=current;
				}else if(last=='D'){
					D++;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m+=current;
		}else if(last=='S'){
			S+=current;
		}else if(last=='C'){
			C+=current;
		}else if(last=='N'){
			N+=current;
		}else if(last=='I'){
			I+=current;
		}else if(last=='D'){
			D++;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		int errors=S+I+D;
		//		if(errors>30) {System.err.println(errors+": "+new String(match));}
		return errors;
	}

	/** Counts m/S/C/N/I/D lengths from expanded symbols or numeric runs.
	 * Requires nonnull input and asserts this six-symbol alphabet.
	 * @return A new array in m,S,C,N,I,D order */
	public static int[] countMatchSymbols(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m+=current;
				}else if(last=='S'){
					S+=current;
				}else if(last=='C'){
					C+=current;
				}else if(last=='N'){
					N+=current;
				}else if(last=='I'){
					I+=current;
				}else if(last=='D'){
					D+=current;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m+=current;
		}else if(last=='S'){
			S+=current;
		}else if(last=='C'){
			C+=current;
		}else if(last=='N'){
			N+=current;
		}else if(last=='I'){
			I+=current;
		}else if(last=='D'){
			D+=current;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return new int[]{m, S, C, N, I, D};
	}

	/**
	 * Counts encoded symbol tokens, ignoring numeric run lengths.
	 * Does not collapse adjacent equal symbols: mmmDDDmmmm yields 7 m and 3 D,
	 * whereas m3D3m4 yields 2 m and 1 D. Compress first for contiguous-run event counts.
	 * @return {m,S,C,N,I,D};
	 */
	public static int[] countMatchEvents(byte[] match){
		int m=0, S=0, C=0, N=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				current=Tools.max(current, 1);
				if(last=='m'){
					m++;
				}else if(last=='S'){
					S++;
				}else if(last=='C'){
					C++;
				}else if(last=='N'){
					N++;
				}else if(last=='I'){
					I++;
				}else if(last=='D'){
					D++;
				}else{
					assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
				}
				current=0;
				last=b;
			}
		}
		current=Tools.max(current, 1);
		if(last=='m'){
			m++;
		}else if(last=='S'){
			S++;
		}else if(last=='C'){
			C++;
		}else if(last=='N'){
			N++;
		}else if(last=='I'){
			I++;
		}else if(last=='D'){
			D++;
		}else{
			assert(last=='?') : "Unhandled symbol "+(char)last+"\n"+new String(match);
		}
		current=0;
		return new int[]{m, S, C, N, I, D};
	}

	/** Sums lengths of uppercase S tokens; requires nonnull input, ignoring other modes. */
	public static int countSubs(byte[] match){
		int S=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				if(last=='S'){S+=Tools.max(1, current);}
				current=0;
				last=b;
			}
		}
		if(last=='S'){S+=Tools.max(1, current);}
		//		assert(S==0) : S+"\t"+new String(match);
		return S;
		//		int x=0;
		//		assert(match!=null);
		//		for(int i=0; i<match.length; i++){
		//			byte b=match[i];
		//			if(b=='S'){x++;}
		//			assert(!Tools.isDigit(b));
		//		}
		//		return x;
	}

	/** Counts all enabled substitution lengths and encoded insertion/deletion tokens. */
	public static int countVars(byte[] match){return countVars(match, true, true, true);}

	/** Sums S lengths and one per I/D token according to the three switches.
	 * Numeric I/D lengths do not multiply event counts; adjacent equal letters are
	 * separate tokens, so compress first for contiguous indel events. Requires nonnull match. */
	public static int countVars(byte[] match, boolean sub, boolean ins, boolean del){
		int S=0, I=0, D=0;
		int current=0;
		byte last='?';
		for(byte b : match){
			if(Tools.isDigit(b)){
				current=current*10+b-'0';
			}else{
				if(last=='S'){S+=Tools.max(1, current);}else if(last=='I'){I++;}else if(last=='D'){D++;}
				current=0;
				last=b;
			}
		}
		if(last=='S'){S+=Tools.max(1, current);}else if(last=='I'){I++;}else if(last=='D'){D++;}
		return (sub ? S : 0)+(ins ? I : 0)+(del ? D : 0);
	}

	/** Tests for a literal uppercase S byte; requires nonnull input. */
	public static boolean containsSubs(byte[] match){
		int x=0;
		assert(match!=null);
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='S'){return true;}
		}
		return false;
	}

	/** Tests for literal S/I/D bytes; requires nonnull input. */
	public static boolean containsVars(byte[] match){
		int x=0;
		assert(match!=null);
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='S' || b=='I' || b=='D'){return true;}
		}
		return false;
	}

	/** Counts literal uppercase N bytes in nonnull input; does not interpret encoded run lengths. */
	public static int countNocalls(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='N'){n++;}
		}
		return n;
	}

	/** Counts nucleotide-undefined bytes over the full nonnull array, independently of any Read amino flag. */
	public static int countUndefined(byte[] bases){return countUndefined(bases, 0, bases.length-1);}

	/** Counts nucleotide-undefined bytes in the inclusive range from..to.
	 * Null input or from&gt;to returns zero; otherwise both endpoints must index bases. */
	public static int countUndefined(byte[] bases, int from, int to){
		int n=0;
		if(bases==null){return 0;}
		for(int i=from; i<=to; i++){
			n+=(AminoAcid.isFullyDefined(bases[i]) ? 0 : 1);
		}
		return n;
	}

	/** Applies the ranged nucleotide homopolymer helper to the full nonnull array. */
	public static int longestHomopolymer(byte[] bases){return longestHomopolymer(bases, 0, bases.length-1);}

	/** Returns the longest exact-byte run of defined nucleotides in inclusive from..to.
	 * Undefined symbols reset the streak to one, so a nonempty all-undefined range returns one.
	 * Case variants do not share a run. Null input or from&gt;to returns zero; otherwise endpoints
	 * must index bases. This helper always uses nucleotide definedness. */
	public static int longestHomopolymer(byte[] bases, int from, int to){
		int max=0, streak=0;
		if(bases==null){return 0;}
		for(int i=from, prev=-1; i<=to; i++){
			final byte b=bases[i];
			if(b==prev && AminoAcid.isFullyDefined(b)){
				streak++;
			}else{
				max=Tools.max(max, streak);
				streak=1;
				prev=b;
			}
		}
		return Tools.max(max, streak);
	}

	/** Counts literal I bytes, not encoded insertion lengths; requires nonnull input. */
	public static int countInsertions(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I'){n++;}
		}
		return n;
	}

	/** Counts literal D bytes, not encoded deletion lengths; requires nonnull input. */
	public static int countDeletions(byte[] match){
		int n=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='D'){n++;}
		}
		return n;
	}

	/** Counts transitions into literal I runs; digit bytes also break runs. */
	public static int countInsertionEvents(byte[] match){
		int n=0;
		byte prev='N';
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='I' && prev!=b){n++;}
			prev=b;
		}
		return n;
	}

	/** Counts transitions into literal D runs; digit bytes also break runs. */
	public static int countDeletionEvents(byte[] match){
		int n=0;
		byte prev='N';
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='D' && prev!=b){n++;}
			prev=b;
		}
		return n;
	}

	/** Tests for literal X/Y anywhere; null returns false. */
	public static boolean containsXY(byte[] match){
		if(match==null){return false;}
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(b=='X' || b=='Y'){return true;}
		}
		return false;
	}

	/** Multiplies PROB_CORRECT for defined bases in the selected prefix, using float arithmetic.
	 * Null qualities return zero. Undefined bases return zero when counted, otherwise
	 * contribute no factor. An empty selected range returns one. Requires valid quality
	 * table indexes and matching positions; maxBases below one selects all qualities. */
	public static float probabilityErrorFree(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		final float[] array=QualityTools.PROB_CORRECT;
		assert(array[0]>0 && array[0]<1);
		float product=1;
		for(int i=0; i<limit; i++){
			byte b=bases[i];
			byte q=quality[i];
			if(AminoAcid.isFullyDefined(b)){
				product*=array[q];
			}else if(countUndefined){
				return 0;
			}
		}
		return product;
	}

	/** Sums PROB_ERROR over a selected prefix; null qualities contribute zero.
	 * Undefined bases contribute the Q0 table value when counted, otherwise nothing.
	 * maxBases below one selects all qualities; requires valid table indexes and matching positions. */
	public static float expectedErrors(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
		final float[] array=QualityTools.PROB_ERROR;
		assert(array[0]>0 && array[0]<1);
		float sum=0;
		for(int i=0; i<limit; i++){
			byte b=bases[i];
			boolean d=AminoAcid.isFullyDefined(b);
			//assert((quality[i]==0)==d) : "Q="+quality[i]+" for base "+(char)b;
			if(d || countUndefined){
				byte q=(d ? quality[i] : 0);
				sum+=array[q];
			}
		}
		return sum;
	}

	/** Sums modeled errors from the quality array's right end; null qualities return zero.
	 * Undefined bases assert Q0 and contribute 0.75 only when countUndefined is true.
	 * maxBases below one selects all qualities; requires matching base positions and valid scores. */
	public static float expectedTipErrors(byte[] bases, byte[] quality, boolean countUndefined, int maxBases){
		if(quality==null){return 0;}
		final int limit;
		{
			final int limit0=(maxBases<1 ? quality.length : Tools.min(maxBases, quality.length));
			limit=quality.length-limit0;
		}
		final float[] array=QualityTools.PROB_ERROR;
		assert(array[0]>0 && array[0]<1);
		float sum=0;
		for(int i=quality.length-1; i>=limit; i--){
			byte b=bases[i];
			byte q=quality[i];
			if(AminoAcid.isFullyDefined(b)){
				sum+=array[q];
			}else{
				assert(q==0);
				if(countUndefined){sum+=0.75f;}
			}
		}
		return sum;
	}

	/** Historical representation heuristic: true for any digit or a run of at least five equal bytes.
	 * Thus an expanded repeated-symbol string may also return true; this is not a strict format test. */
	public static boolean isShortMatchString(byte[] match){
		byte last=' ';
		int streak=0;
		for(int i=0; i<match.length; i++){
			byte b=match[i];
			if(Tools.isDigit(b)){return true;}
			if(b==last){
				streak++;
				if(streak>3){return true;}
			}else{
				streak=0;
				last=b;
			}
		}
		return false;
	}

	/** Run-length encodes expanded letter bytes into a new array; null remains null.
	 * Counts greater than one follow their letter. Requires nonempty input under assertions;
	 * does not update any Read's representation flag. */
	public static byte[] toShortMatchString(byte[] match){
		if(match==null){return null;}
		assert(match.length>0);
		ByteBuilder sb=new ByteBuilder(10);

		byte prev=match[0];
		int count=1;
		for(int i=1; i<match.length; i++){
			byte m=match[i];
			assert(Tools.isLetter(m) || m==0) : new String(match);
			if(m==0){System.err.println("Warning! Converting empty match string to short form.");}
			if(m==prev){count++;}else{
				sb.append(prev);
				if(count>1){sb.append(count);}
				//				else if(count==2){sb.append(prev);}
				prev=m;
				count=1;
			}
		}
		sb.append(prev);
		if(count>1){sb.append(count);}
		//		else if(count==2){sb.append(prev);}

		byte[] r=sb.toBytes();
		return r;
	}

	//STUDIED-PRAISE (with toShortMatchString): the RLE match codec. This decoder runs TWO passes - the first computes the EXACT decompressed length so the second fills one right-sized array with no resizing. The per-letter length formula `count += (current>0 ? current-1 : 0)` exactly mirrors the fill loop's `while(current>1)`, so the passes cannot disagree (a mismatch would AIOOBE loudly, never silently corrupt). Round-trip verified: toShort("mmmmSS")="m4S2", toLong back to "mmmmSS".
	/** Decodes letter/count runs into a newly allocated expanded array, preserving null.
	 * The caller supplies a nonempty letter/count sequence; assertions check only some
	 * constraints. Does not update a Read flag.
	 * Two passes size and fill the result using the same count convention. */
	public static byte[] toLongMatchString(byte[] shortmatch){
		if(shortmatch==null){return null;}
		assert(shortmatch.length>0);

		int count=0;
		int current=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Tools.isLetter(m)){
				count++;
				count+=(current>0 ? current-1 : 0);
				current=0;
			}else{
				assert(Tools.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		count+=(current>0 ? current-1 : 0);


		byte[] r=new byte[count];
		current=0;
		byte lastLetter='?';
		int j=0;
		for(int i=0; i<shortmatch.length; i++){
			byte m=shortmatch[i];
			if(Tools.isLetter(m)){
				while(current>1){
					r[j]=lastLetter;
					current--;
					j++;
				}
				current=0;

				r[j]=m;
				j++;
				lastLetter=m;
			}else{
				assert(Tools.isDigit(m));
				current=(current*10)+(m-48); //48 == '0'
			}
		}
		while(current>1){
			r[j]=lastLetter;
			current--;
			j++;
		}

		assert(r[r.length-1]>0);
		return r;
	}

	/** Selects the plus-left estimate only for two mapped opposite-strand reads when requested.
	 * Otherwise uses the unstranded estimate, including when a mapped bit is clear.
	 * Requires nonnull r1; r2 may be null. Does not validate alignment coordinates. */
	public static int insertSizeMapped(Read r1, Read r2, boolean ignoreStrand){
		//		assert(false) : ignoreStrand+", "+(r2==null)+", "+(r1.mapped())+", "+(r2.mapped())+", "+(r1.strand()==r2.strand())+", "+r1.strand()+", "+r2.strand();
		if(ignoreStrand || r2==null || !r1.mapped() || !r2.mapped() || r1.strand()==r2.strand()){return insertSizeMapped_Unstranded(r1, r2);}
		return insertSizeMapped_PlusLeft(r1, r2);
	}

	/** Orders two nonnull reads by strand and estimates gap plus both sequence lengths.
	 * Falls back to the unstranded estimate for same strands, reversed endpoints or
	 * excessive overlap; different chromosomes or a singleton coordinate span return zero.
	 * Does not test mapped flags. TODO: This is not correct when the insert is shorter
	 * than a read's bases with same-strand reads. */
	public static int insertSizeMapped_PlusLeft(Read r1, Read r2){
		if(r1.strand()>r2.strand()){return insertSizeMapped_PlusLeft(r2, r1);}
		if(r1.strand()==r2.strand() || r1.start>r2.stop){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left.
		//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.chrom!=r2.chrom){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???

		int a=r1.length();
		int b=r2.length();
		int mid=r2.start-r1.stop-1;
		if(-mid>=a+b){return insertSizeMapped_Unstranded(r1, r2);} //Not properly oriented; plus read is to the right of minus read
		return mid+a+b;
	}

	/** Estimates insert size without strand or mapped-flag checks; requires nonnull r1.
	 * A missing r2 yields r1's inclusive coordinate span, except singleton spans yield zero.
	 * With two reads, different chromosomes or singleton spans yield zero; equal starts
	 * yield the shorter sequence length. Other cases use gap plus both lengths, returning
	 * zero when the overlap is at least their summed sequence length. */
	public static int insertSizeMapped_Unstranded(Read r1, Read r2){
		if(r2==null){return r1.start==r1.stop ? 0 : r1.stop-r1.start+1;}

		if(r1.start>r2.start){return insertSizeMapped_Unstranded(r2, r1);} //So r1 is always on the left side.

		//		if(!mapped() || !mate.mapped()){return 0;}
		if(r1.start==r1.stop || r2.start==r2.stop){return 0;} //???

		if(r1.chrom!=r2.chrom){return 0;}
		int a=r1.length();
		int b=r2.length();
		if(false && Tools.overlap(r1.start, r1.stop, r2.start, r2.stop)){
			//This does not handle very short inserts
			return Tools.max(r1.stop, r2.stop)-Tools.min(r1.start, r2.start)+1;

		}else{
			if(r1.start<r2.start){
				int mid=r2.start-r1.stop-1;
				//				assert(false) : mid+", "+a+", "+b;
				//				if(-mid>a && -mid>b){return Tools.min(a, b);} //Strange situation, no way to guess insert size
				if(-mid>=a+b){return 0;} //Strange situation, no way to guess insert size
				return mid+a+b;
			}else{
				assert(r1.start==r2.start);
				return Tools.min(a, b);
			}
		}
	}

	/** Applies the coordinate insert estimate to two nonnull sites and supplied read lengths. */
	public static int insertSize(SiteScore ssa, SiteScore ssb, int lena, int lenb){return insertSize(ssa.chrom, ssb.chrom, ssa.start, ssb.start, ssa.stop, ssb.stop, lena, lenb);}

	/** Returns zero for different chromosomes, the outer span for overlapping intervals,
	 * or their intervening gap plus lena+lenb otherwise. Endpoints are inclusive;
	 * strands, mapped flags and consistency between lengths and coordinates are not checked. */
	public static int insertSize(int chroma, int chromb, int starta, int startb, int stopa, int stopb, int lena, int lenb){

		final int x;

		//		if(mate==null || ){return bases==null ? 0 : bases.length;}
		if(chroma!=chromb){x=0;}else{

			if(Tools.overlap(starta, stopa, startb, stopb)){
				x=Tools.max(stopa, stopb)-Tools.min(starta, startb)+1;
				//				System.err.println("C: "+x);
			}else{
				if(starta<=startb){
					int mid=startb-stopa-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
					//					System.err.println("D: "+x);
				}else{
					int mid=starta-stopb-1;
					//				assert(false) : mid+", "+a+", "+b;
					x=mid+lena+lenb;
					//					System.err.println("E: "+x);
				}
			}
		}
		return x;
	}

	/** Creates an insert-length read from already-oriented, nonnull reads with nonnull bases.
	 * Does not reverse-complement either input. Copies a from the left and b from the right,
	 * filling a gap with N/Q0. Overlap with qualities favors the higher-quality call,
	 * uses N on equal-quality disagreements, and boosts agreement quality up to the merge cap.
	 * Without qualities, replaces empty/N positions from b and otherwise favors the larger byte.
	 * Output qualities are absent unless both inputs supply them; see the asymmetric assertion below.
	 * The constructor may validate the new bases before joined qualities are assigned.
	 * Starts with a's ID/flags/chromosome; the constructor clears its validation bit before optional validation.
	 * Then sets insert and clears the paired bit; no mate is linked.
	 * Can set a.chrom to zero when clearing output coordinates. Legacy mapped-bit behavior and
	 * quality-presence assumptions remain separately annotated, not repaired here.
	 * @param a Left-oriented input, whose chromosome may be cleared
	 * @param b Right-oriented input
	 * @param insert Positive output length
	 * @return New read with independently allocated base/optional quality arrays */
	public static Read joinRead(Read a, Read b, int insert){
		assert(a!=null && b!=null && insert>0);
		final int lengthSum=a.length()+b.length();
		final int overlap=Tools.min(insert, lengthSum-insert);

		//		System.err.println(insert);
		final byte[] bases=new byte[insert], abases=a.bases, bbases=b.bases;
		final byte[] aquals=a.quality, bquals=b.quality;
		final byte[] quals=(aquals==null || bquals==null ? null : new byte[insert]);
		//TODO: Quality-presence asymmetry - nonnull, length-matched aquals dereferences bquals here,
		//although allocation above accepts either side absent. Mixed-quality caller support is unverified.
		assert(aquals==null || (aquals.length==abases.length && bquals.length==bbases.length));

		int mismatches=0;

		int start, stop;

		if(overlap<=0){//Simple join in which there is no overlap
			int lim=insert-b.length();
			if(quals==null){
				for(int i=0; i<a.length(); i++){
					bases[i]=abases[i];
				}
				for(int i=a.length(); i<lim; i++){
					bases[i]='N';
				}
				for(int i=0; i<b.length(); i++){
					bases[i+lim]=bbases[i];
				}
			}else{
				for(int i=0; i<a.length(); i++){
					bases[i]=abases[i];
					quals[i]=aquals[i];
				}
				for(int i=a.length(); i<lim; i++){
					bases[i]='N';
					quals[i]=0;
				}
				for(int i=0; i<b.length(); i++){
					bases[i+lim]=bbases[i];
					quals[i+lim]=bquals[i];
				}
			}

			start=Tools.min(a.start, b.start);
			//			stop=start+insert-1;
			stop=Tools.max(a.stop, b.stop);

			//		}else if(insert>=a.length() && insert>=b.length()){ //Overlapped join, proper orientation
			//			final int lim1=a.length()-overlap;
			//			final int lim2=a.length();
			//			for(int i=0; i<lim1; i++){
			//				bases[i]=abases[i];
			//				quals[i]=aquals[i];
			//			}
			//			for(int i=lim1, j=0; i<lim2; i++, j++){
			//				assert(false) : "TODO";
			//				bases[i]='N';
			//				quals[i]=0;
			//			}
			//			for(int i=lim2, j=overlap; i<bases.length; i++, j++){
			//				bases[i]=bbases[j];
			//				quals[i]=bquals[j];
			//			}
		}else{ //reads go off ends of molecule.
			if(quals==null){
				for(int i=0; i<a.length() && i<bases.length; i++){
					bases[i]=abases[i];
				}
				for(int i=bases.length-1, j=b.length()-1; i>=0 && j>=0; i--, j--){
					byte ca=bases[i], cb=bbases[j];
					if(ca==0 || ca=='N'){
						bases[i]=cb;
					}else if(ca==cb){
					}else{
						bases[i]=(ca>=cb ? ca : cb);
						if(ca!='N' && cb!='N'){mismatches++;}
					}
				}
			}else{
				for(int i=0; i<a.length() && i<bases.length; i++){
					bases[i]=abases[i];
					quals[i]=aquals[i];
				}
				for(int i=bases.length-1, j=b.length()-1; i>=0 && j>=0; i--, j--){
					byte ca=bases[i], cb=bbases[j];
					byte qa=quals[i], qb=bquals[j];
					if(ca==0 || ca=='N'){
						bases[i]=cb;
						quals[i]=qb;
					}else if(cb==0 || cb=='N'){
						//do nothing
					}else if(ca==cb){
						quals[i]=(byte)Tools.min((Tools.max(qa, qb)+Tools.min(qa, qb)/4), MAX_MERGE_QUALITY);
					}else{
						bases[i]=(qa>qb ? ca : qa<qb ? cb : (byte)'N');
						quals[i]=(byte)(Tools.max(qa, qb)-Tools.min(qa, qb));
						if(ca!='N' && cb!='N'){mismatches++;}
					}
				}
			}

			if(a.strand()==0){
				start=a.start;
				//				stop=start+insert-1;
				stop=b.stop;
			}else{
				stop=a.stop;
				//				start=stop-insert+1;
				start=b.start;
			}
			if(start>stop){
				start=Tools.min(a.start, b.start);
				stop=Tools.max(a.stop, b.stop);
			}
		}
		//		assert(mismatches>=countMismatches(a, b, insert, 999));
		//		System.err.println(mismatches);
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){start=stop=a.chrom=0;}

		//		System.err.println(bases.length+", "+start+", "+stop);

		Read r=new Read(bases, null, a.id, a.numericID, a.flags, a.chrom, start, stop);
		r.quality=quals; //This prevents quality from getting capped.
		//TODO: Possible bug [stream/Read#006] - QUESTION: sets r.setMapped(TRUE) under the SAME condition that just zeroed the coords (~L3859, the "unmapped-ish" case). Marking a coord-zeroed read MAPPED looks inverted (setMapped(false)?). Hot BBMerge code - NOT auto-fixed; FLAG FOR BRIAN. Possible MEDIUM (wrong mapped flag on merged output) if unintentional.
		if(a.chrom==0 || start==stop || (!a.mapped() && !a.synthetic())){r.setMapped(true);}
		r.setInsert(insert);
		r.setPaired(false);
		r.copies=a.copies;
		r.mapScore=a.mapScore+b.mapScore;
		if(overlap<=0){
			r.mapScore=a.mapScore+b.mapScore;
			r.errors=a.errors+b.errors;
			//TODO r.gaps=?
		}else{//Hard to calculate
			r.mapScore=(int)((a.mapScore*(long)a.length()+b.mapScore*(long)b.length())/insert);
			r.errors=a.errors;
		}


		assert(r.insertvalid()) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		assert(r.insert()==r.length()) : r.insert()+"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";
		//		assert(false) : "\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";

		//TODO: Triggered by BBMerge in useratio mode for some reason.
		//		assert(Shared.anomaly || (a.insertSizeMapped(false)>0 == r.insertSizeMapped(false)>0)) :
		//			"\n"+r.length()+"\n"+r.insert()+"\n"+r.insertSizeMapped(false)+"\n"+a.insert()+"\n"+a.insertSizeMapped(false)+
		//			"\n\n"+a.toText(false)+"\n\n"+b.toText(false)+"\n\n"+r.toText(false)+"\n\n";

		return r;
	}

	//	/** Generate and return an array of canonical kmers for this read */
	//	public long[] toLongKmers(final int k, long[] kmers, boolean makeCanonical, Kmer longkmer) {
	//		assert(k>31) : k;
	//		if(bases==null || bases.length<k){return null;}
	//
	//		final int kbits=2*k;
	//		final long mask=Long.MAX_VALUE;
	//
	//		int len=0;
	//		long kmer=0;
	//		final int arraylen=bases.length-k+1;
	//		if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
	//		Arrays.fill(kmers, -1);
	//
	//
	//		final int tailshift=k%32;
	//		final int tailshiftbits=tailshift*2;
	//
	//		for(int i=0; i<bases.length; i++){
	//			byte b=bases[i];
	//			int x=AminoAcid.baseToNumber[b];
	//			if(x<0){
	//				len=0;
	//				kmer=0;
	//			}else{
	//				kmer=Long.rotateLeft(kmer, 2);
	//				kmer=kmer^x;
	//				len++;
	//
	//				if(len>=k){
	//					long x2=AminoAcid.baseToNumber[bases[i-k]];
	//					kmer=kmer^(x2<<tailshiftbits);
	//					kmers[i-k+1]=kmer;
	//				}
	//			}
	//		}
	//		if(makeCanonical){
	//			this.reverseComplement();
	//			len=0;
	//			kmer=0;
	//			for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
	//				byte b=bases[i];
	//				int x=AminoAcid.baseToNumber[b];
	//				if(x<0){
	//					len=0;
	//					kmer=0;
	//				}else{
	//					kmer=Long.rotateLeft(kmer, 2);
	//					kmer=kmer^x;
	//					len++;
	//
	//					if(len>=k){
	//						long x2=AminoAcid.baseToNumber[bases[i-k]];
	//						kmer=kmer^(x2<<tailshiftbits);
	//						kmers[j]=mask&(Tools.max(kmers[j], kmer));
	//					}
	//				}
	//			}
	//			this.reverseComplement();
	//		}else{
	//			assert(false) : "Long kmers should be made canonical here because they cannot be canonicized later.";
	//		}
	//
	//		return kmers;
	//	}

	/** Forwards this read's fields to the currently disabled site checker. */
	public static final boolean CHECKSITES(Read r, byte[] basesM){return CHECKSITES(r.sites, r.bases, basesM, r.numericID, true);}

	/** Forwards this read's fields and requested sorting check to the disabled checker. */
	public static final boolean CHECKSITES(Read r, byte[] basesM, boolean verifySorted){return CHECKSITES(r.sites, r.bases, basesM, r.numericID, verifySorted);}

	/** Delegates to the disabled checker with verifySorted=true; performs no validation. */
	public static final boolean CHECKSITES(ArrayList<SiteScore> list, byte[] basesP, byte[] basesM, long id){return CHECKSITES(list, basesP, basesM, id, true);}

	/** Disabled historical checker: always returns true and ignores the supplied state. */
	public static final boolean CHECKSITES(ArrayList<SiteScore> list, byte[] basesP, byte[] basesM, long id, boolean verifySorted){
		return true; //Temporarily disabled
		//		if(list==null || list.isEmpty()){return true;}
		//		SiteScore prev=null;
		//		for(int i=0; i<list.size(); i++){
		//			SiteScore ss=list.get(i);
		//			if(ss.strand==Gene.MINUS && basesM==null && basesP!=null){basesM=AminoAcid.reverseComplementBases(basesP);}
		//			byte[] bases=(ss.strand==Gene.PLUS ? basesP : basesM);
		//			if(verbose){System.err.println("Checking site "+i+": "+ss);}
		//			boolean b=CHECKSITE(ss, bases, id);
		//			assert(b) : id+"\n"+new String(basesP)+"\n"+ss+"\n";
		//			if(verbose){System.err.println("Checked site "+i+" = "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect);}
		//			if(!b){
		////				System.err.println("Error at SiteScore "+i+": ss.p="+ss.perfect+", ss.sp="+ss.semiperfect);
		//				return false;
		//			}
		//			if(verifySorted && prev!=null && ss.score>prev.score){
		//				if(verbose){System.err.println("verifySorted failed.");}
		//				return false;
		//			}
		//			prev=ss;
		//		}
		//		return true;
	}

	/** Tests nonincreasing score order only; null and lists shorter than two pass. */
	public static final boolean CHECKORDER(ArrayList<SiteScore> list){
		if(list==null || list.size()<2){return true;}
		SiteScore prev=list.get(0);
		for(int i=0; i<list.size(); i++){
			SiteScore ss=list.get(i);
			if(ss.score>prev.score){return false;}
			prev=ss;
		}
		return true;
	}

	/** Selects the strand-oriented base array, then delegates to the disabled checker. */
	public static final boolean CHECKSITE(SiteScore ss, byte[] basesP, byte[] basesM, long id){return CHECKSITE(ss, ss.plus() ? basesP : basesM, id);}

	/** Disabled historical checker: returns true without inspecting site, bases or ID. */
	public static final boolean CHECKSITE(SiteScore ss, byte[] bases, long id){
		return true; //Temporarily disabled
		//		if(ss==null){return true;}
		////		System.err.println("Checking site "+ss+"\nss.p="+ss.perfect+", ss.sp="+ss.semiperfect+", bases="+new String(bases));
		//		if(ss.perfect){assert(ss.semiperfect) : ss+"\n"+new String(bases);}
		//		if(ss.gaps!=null){
		//			if(ss.gaps[0]!=ss.start || ss.gaps[ss.gaps.length-1]!=ss.stop){return false;}
		////			assert(ss.gaps[0]==ss.start && ss.gaps[ss.gaps.length-1]==ss.stop);
		//		}
		//
		//		if(!(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore)){
		//			System.err.println("Site paired score violation: "+ss.quickScore+", "+ss.slowScore+", "+ss.pairedScore);
		//			return false;
		//		}
		//
		//		final boolean xy=ss.matchContainsXY();
		//		if(bases!=null){
		//
		//			final boolean p0=ss.perfect;
		//			final boolean sp0=ss.semiperfect;
		//			final boolean p1=ss.isPerfect(bases);
		//			final boolean sp1=(p1 ? true : ss.isSemiPerfect(bases));
		//
		//			assert(p0==p1 || (xy && p1)) : p0+"->"+p1+", "+sp0+"->"+sp1+", "+ss.isSemiPerfect(bases)+
		//				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n"+ss+"\n\n";
		//			assert(sp0==sp1 || (xy && sp1)) : p0+"->"+p1+", "+sp0+"->"+sp1+", "+ss.isSemiPerfect(bases)+
		//				"\nnumericID="+id+"\n"+new String(bases)+"\n\n"+Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n\n"+ss+"\n\n";
		//
		////			ss.setPerfect(bases, false);
		//
		//			assert(p0==ss.perfect) :
		//				p0+"->"+ss.perfect+", "+sp0+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n\n"+new String(bases)+"\n\n"+
		//				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n"+ss+"\n\n";
		//			assert(sp0==ss.semiperfect) :
		//				p0+"->"+ss.perfect+", "+sp0+"->"+ss.semiperfect+", "+ss.isSemiPerfect(bases)+"\nnumericID="+id+"\n\n"+new String(bases)+"\n\n"+
		//				Data.getChromosome(ss.chrom).getString(ss.start, ss.stop)+"\n"+ss+"\n\n";
		//			if(ss.perfect){assert(ss.semiperfect);}
		//		}
		//		if(ss.match!=null && ss.matchLength()!=ss.mappedLength()){
		//			if(verbose){System.err.println("Returning false because matchLength!=mappedLength:\n"+ss.matchLength()+", "+ss.mappedLength()+"\n"+ss);}
		//			return false;
		//		}
		//		return true;
	}

	/** Allocates max+1 masks with bit i set at index i; Java int shift-distance rules apply. */
	private static int[] makeMaskArray(int max){
		int[] r=new int[max+1];
		for(int i=0; i<r.length; i++){r[i]=(1<<i);}
		return r;
	}



	/** Returns Q30 bytes: a shared mutable cache entry below length 1000, otherwise a new array.
	 * Treat cached results as borrowed; this is not an independent quality buffer.
	 * Cache publication has the separate source concern below. */
	public static byte[] getFakeQuality(int len){
		if(len>=QUALCACHE.length){
			byte[] r=KillSwitch.allocByte1D(len);
			Arrays.fill(r, (byte)30);
			return r;
		}
		if(QUALCACHE[len]==null){
			synchronized(QUALCACHE){
				if(QUALCACHE[len]==null){
					//TODO: Probable concurrency bug - publish-before-fill lets an unsynchronized
					//reader above obtain this array before Q30 initialization is complete.
					QUALCACHE[len]=KillSwitch.allocByte1D(len);
					Arrays.fill(QUALCACHE[len], (byte)30);
				}
			}
		}
		return QUALCACHE[len];
	}

	/** Returns the current global minimum called-quality setting. */
	public static byte MIN_CALLED_QUALITY(){return MIN_CALLED_QUALITY;}
	/** Returns the current global maximum called-quality setting, independently of qMap contents. */
	public static byte MAX_CALLED_QUALITY(){return MAX_CALLED_QUALITY;}

	/** Restricts x to 1..93 and refills the shared lookup only if the maximum changes.
	 * Does not adjust the minimum; this is global mutable configuration, not per-read state. */
	public static void setMaxCalledQuality(int x){
		x=Tools.mid(1, x, 93);
		if(x!=MAX_CALLED_QUALITY){
			MAX_CALLED_QUALITY=(byte)x;
			qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);
		}
	}

	/** Restricts x to 0..93 and refills the shared lookup only if the minimum changes.
	 * Does not adjust the maximum; callers must coordinate use of the mutable lookup. */
	public static void setMinCalledQuality(int x){
		x=Tools.mid(0, x, 93);
		if(x!=MIN_CALLED_QUALITY){
			MIN_CALLED_QUALITY=(byte)x;
			qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);
		}
	}

	/** Returns the median of q and the current global minimum/maximum, narrowed to byte. */
	public static byte capQuality(long q){return (byte)Tools.mid(MIN_CALLED_QUALITY, q, MAX_CALLED_QUALITY);}

	/** Looks up q in the shared table; requires a valid nonnegative index and performs no clamping first. */
	public static byte capQuality(byte q){return qMap[q];}

	/** Looks up q for a fully defined base, otherwise returns zero without accessing qMap.
	 * Defined bases require q to be a valid table index. */
	public static byte capQuality(byte q, byte b){return AminoAcid.isFullyDefined(b) ? qMap[q] : 0;}

	/** Fills qMap in place, or allocates 128 entries if null, with median(min,index,max).
	 * Returns the same mutable array it fills; supplied limits are not reordered or stored here. */
	private static byte[] makeQmap(byte min, byte max){
		byte[] r=(qMap==null ? new byte[128] : qMap);
		for(int i=0; i<r.length; i++){
			r[i]=(byte) Tools.mid(min, i, max);
		}
		return r;
	}

	/*--------------------------------------------------------------*/
	/*--------------------        Fields        --------------------*/
	/*--------------------------------------------------------------*/

	/** Serialization version retained for compatibility with saved Read objects. */
	private static final long serialVersionUID=-1026645233407290096L;

	/** Mutable base/AA symbols; constructors retain the caller's array. */
	public byte[] bases;

	/** Mutable numeric Phred qualities, or null; not ASCII-encoded FASTQ characters. */
	public byte[] quality;

	/** Alignment string.  E.G. mmmmDDDmmm would have 4 matching bases, then a 3-base deletion, then 3 matching bases. */
	public byte[] match;

	/** Gap coordinates for alignments spanning multiple reference regions */
	public int[] gaps;

	/** Read identifier string from sequence header */
	public String id;
	/** Numeric identifier for efficient processing */
	public long numericID;
	/** Reference chromosome/scaffold ID for mapped reads */
	public int chrom=-1;
	/** Mapping start position (0-based) */
	public int start=-1;
	/** Mapping stop position (0-based, inclusive) */
	public int stop=-1;

	/** Number of duplicate copies merged into this read */
	public int copies=1;

	/** Errors detected */
	public int errors=0;

	/** Alignment score from BBMap.  Assumed to max at approx 100*bases.length */
	public int mapScore=0;
	/** Runtime-only packed loose/strict neural MAPQs; zero means neither is present. */
	public short neuralMapqs=0;

	/** Mutable candidate sites, conventionally score-sorted by callers; not sorted on assignment. */
	public ArrayList<SiteScore> sites;
	/** Original alignment site for synthetic reads */
	public SiteScore originalSite; //Origin site for synthetic reads
	/** Caller-specific payload shared by shallow copies; use samline for a SamLine. */
	public Object obj=null;//Don't set this to a SamLine; it's for other things.
	/** Retained SAM representation; field mutations do not automatically regenerate it. */
	public SamLine samline=null;
	/** Paired-end mate read */
	public Read mate;

	/**
	 * Bit flags storing read properties (strand, mapping status, quality flags, etc.)
	 */
	public int flags;

	/** -1 if invalid.  TODO: Currently not retained through most processes. */
	private int insert=-1;

	/** A random number for deterministic usage.
	 * May decrease speed in multithreaded applications.
	 */
	public double rand=-1;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields        -----------------*/
	/*--------------------------------------------------------------*/

	/** Shared lazily filled Q30 arrays indexed by length; callers must treat entries as borrowed. */
	private static final byte[][] QUALCACHE=new byte[1000][];


	//Flag bit-layout (the `flags` int): bit 0=strand, 1=mapped, 2=paired, 3=perfect, 4=ambiguous, 5=rescued, 6=supplementary (old COLORMASK), 7=synthetic, 8=discard, 9=invalid, 10=swap, 11=shortmatch, 12=pairnum (single bit: 0=read1/1=read2), 13=insert-valid, 14=adapter, 15=secondary, 16=aminoacid, 17=junk, 18=validated, 19=tested, 20=inverted-repeat, 21=trimmed. 22 bits of 32 used; accessors are (flags&MASK)==MASK, setters clear-then-set. toText writes them as a binary string over maskArray (bits 0-22), fromText does Integer.parseInt(...,2) - symmetric.
	/** Low strand bit: zero is plus, one is minus. */
	public static final int STRANDMASK=1;
	/** Stored mapped-state bit; independent of coordinates and mate state. */
	public static final int MAPPEDMASK=(1<<1);
	/** Paired-alignment bit; a mate reference alone does not set it. */
	public static final int PAIREDMASK=(1<<2);
	/** Stored perfect-alignment bit. */
	public static final int PERFECTMASK=(1<<3);
	/** Stored ambiguous-mapping bit. */
	public static final int AMBIMASK=(1<<4);
	/** Stored mate-rescue bit. */
	public static final int RESCUEDMASK=(1<<5);
	/** Supplementary-alignment bit in Read flags; distinct from the SAM bit position. */
	public static final int SUPPLEMENTARYMASK=(1<<6); //SAM flag 0x800; reuses the freed old COLORMASK bit
	/** Synthetic-read bit. */
	public static final int SYNTHMASK=(1<<7);
	/** Discard-request bit; consumers decide whether to omit the read. */
	public static final int DISCARDMASK=(1<<8);
	/** Invalid-state bit tested by valid() and invalid(). */
	public static final int INVALIDMASK=(1<<9);
	/** Stored swap-state bit. */
	public static final int SWAPMASK=(1<<10);
	/** Match-representation bit: set for run-counted short match strings. */
	public static final int SHORTMATCHMASK=(1<<11);

	/** Position of the one-bit read1/read2 selector. */
	public static final int PAIRNUMSHIFT=12;
	/** Mate-number bit: zero identifies read1, one identifies read2. */
	public static final int PAIRNUMMASK=(1<<PAIRNUMSHIFT);

	/** Stored insert-valid bit; separate from the insert value. */
	public static final int INSERTMASK=(1<<13);
	/** Stored adapter-detection bit. */
	public static final int ADAPTERMASK=(1<<14);
	/** Secondary-alignment bit. */
	public static final int SECONDARYMASK=(1<<15);
	/** Amino-acid alphabet bit used by validation and quality normalization. */
	public static final int AAMASK=(1<<16);
	/** Junk-classification bit; separate from invalid, discarded and validated. */
	public static final int JUNKMASK=(1<<17);
	/** Records that validation ran, not that the read passed every check. */
	public static final int VALIDATEDMASK=(1<<18);
	/** Caller-managed tested-state bit. */
	public static final int TESTEDMASK=(1<<19);
	/** Stored inverted-repeat bit. */
	public static final int IRMASK=(1<<20);
	/** Stored trimmed-state bit; not an assertion about the auxiliary payload. */
	public static final int TRIMMEDMASK=(1<<21);

	/** Ascending single-bit masks through bit 22 for native-text flag serialization. */
	private static final int[] maskArray=makeMaskArray(22); //Be sure this is big enough for all flags!

	/** Enables case normalization. fixCase favors this over LOWER_CASE_TO_N; fixCaseVec does the reverse. */
	public static boolean TO_UPPER_CASE=false;
	/** Enables lowercase-to-nocall conversion using the active alphabet and helper precedence. */
	public static boolean LOWER_CASE_TO_N=false;
	/** Enables the active alphabet helper's dot/dash/X-to-nocall conversion. */
	public static boolean DOT_DASH_X_TO_N=false;
	/** Selects probability-derived Phred averages; false selects arithmetic score averages. */
	public static boolean AVERAGE_QUALITY_BY_PROBABILITY=true;
	/** Calls Tools.fixHeader during validate(), including the null-base path. */
	public static boolean FIX_HEADER=false;
	/** Lets header repair replace a null or empty identifier with an empty string. */
	public static boolean ALLOW_NULL_HEADER=false;
	/** Skips sequence normalization and junk checks; quality-length and enabled header checks still run. */
	public static boolean SKIP_SLOW_VALIDATION=false;
	/** Compile-time selection of the branchless nucleotide fallback; the legacy scalar branch is not selected. */
	public static final boolean VALIDATE_BRANCHLESS=true;
	/** Selects vector-dispatched nucleotide validation when Shared.SIMD is also enabled. */
	public static boolean VALIDATE_VECTOR=true;

	/** Skips junk rejection while still allowing separately configured conversions. */
	public static final int IGNORE_JUNK=0;
	/** Selects junk marking rather than replacement or a fatal junk diagnostic. */
	public static final int FLAG_JUNK=1;
	/** Selects junk replacement; the alphabet and active helper determine its conversions. */
	public static final int FIX_JUNK=2;
	/** Selects fatal junk diagnostics when the active helper enables them. */
	public static final int CRASH_JUNK=3;
	/** Selects additional nucleotide-to-ACGTN conversion; vector and branchless trigger conditions differ. */
	public static final int FIX_JUNK_AND_IUPAC=4;
	/** Process-wide junk-handling selector; see the active validation helper for exact behavior. */
	public static int JUNK_MODE=CRASH_JUNK;
	/** Requests IUPAC-to-N conversion in the active nucleotide validation paths. */
	public static boolean IUPAC_TO_N=false;

	/** Converts U/u before nucleotide validation unless SKIP_SLOW_VALIDATION bypasses that path. */
	public static boolean U_TO_T=false;
	/** Allows toText to temporarily short-encode match when its caller also requests compression. */
	public static boolean COMPRESS_MATCH_BEFORE_WRITING=true;
	/** Expands a short-flagged match after native-text parsing. */
	public static boolean DECOMPRESS_MATCH_ON_LOAD=true; //Set to false for some applications, like sorting, perhaps

	/** Selects native-text site parsing from index 17 when true, 18 when false; see the decoder layout TODO. */
	public static boolean ADD_BEST_SITE_TO_LIST_FROM_TEXT=true;
	/** Clears length-mismatched qualities and marks junk, unless TOSS_BROKEN_QUALITY takes precedence. */
	public static boolean NULLIFY_BROKEN_QUALITY=false;
	/** Clears length-mismatched qualities and marks discarded/junk; takes precedence over NULLIFY and FLAG. */
	public static boolean TOSS_BROKEN_QUALITY=false;
	/** Marks a quality-length mismatch as junk while retaining qualities, unless TOSS or NULLIFY is selected. */
	public static boolean FLAG_BROKEN_QUALITY=false;
	/** Selects flat identity scoring; false selects square-root deletion weighting. */
	public static boolean FLAT_IDENTITY=true;
	/** Default for constructors without an explicit validation choice; explicit boolean overloads override it. */
	public static boolean VALIDATE_IN_CONSTRUCTOR=true;

	/** Legacy diagnostic switch; Read's current diagnostic uses are commented out. */
	public static boolean verbose=false;

	/*--------------------------------------------------------------*/

	/** Fixed native-text quality offset; independent of FASTQ's configurable offsets. */
	private static final byte ASCII_OFFSET=33;

	/** Enables ordinary validation quality normalization; repair-specific paths may also change qualities. */
	public static boolean CHANGE_QUALITY=true; //Cap all quality values between MIN_CALLED_QUALITY and MAX_CALLED_QUALITY

	/** Minimum allowed quality score for called (ACGT) bases */
	private static byte MIN_CALLED_QUALITY=2;

	/** Maximum allowed quality score for called (ACGT) bases */
	public static byte MAX_CALLED_QUALITY=50; //TODO: Find, replace, and test all instances of 41 (old value).

	/** Maximum allowed quality score for merged reads, which otherwise would normally be very high */
	public static byte MAX_MERGE_QUALITY=50;
	/** Shared mutable quality lookup, normally 128 entries; setters refill the existing array.
	 * Direct changes to MAX_CALLED_QUALITY do not rebuild it. Coordinate configuration with users. */
	public static byte[] qMap=makeQmap(MIN_CALLED_QUALITY, MAX_CALLED_QUALITY);

}
