package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.SamLine;

/**
 * Utility class for analyzing and manipulating variants in reads.
 * Contains static methods for variant detection, modification, and filtering.
 * Processes genomic reads to identify substitutions, insertions, and deletions,
 * and provides functionality for marking known variants and loading VCF data.
 *
 * @author Brian Bushnell
 * @author Isla Winglet
 * @date December 2024
 */
public class AnalyzeVars {

	/**
	 * Fixes variants in a read by changing match string characters to indicate known variants.
	 * Changes 'S' to 'V' for substitutions, 'I' to 'i' for insertions, 'D' to 'd' for deletions.
	 * Performs initial validation and delegates to the main fixVars implementation.
	 *
	 * @param r Read to modify
	 * @param varMap Map of known variants for comparison
	 * @param scafMap Scaffold mapping information for coordinate resolution
	 * @return Number of variants fixed, 0 if read is invalid or unmapped
	 */
	public static int fixVars(Read r, VarMap varMap, ScafMap scafMap){
		if(r==null || r.bases==null || r.match==null || r.samline==null){return 0;}
		final SamLine sl=r.samline;
		if(!sl.mapped()){return 0;}
		return AnalyzeVars.fixVars(r, sl, varMap, scafMap);
	}

	/**
	 * Reverses the effects of fixVars by changing variant indicators back to standard match characters.
	 * Changes 'V' back to 'S', 'i' back to 'I', 'd' back to 'D'.
	 * Used to restore original match string format after variant processing.
	 * @param r Read to modify, null-safe operation
	 */
	public static void unfixVars(Read r){
		if(r==null || r.match==null){return;}
		for(int i=0; i<r.match.length; i++){
			if(r.match[i]=='V'){r.match[i]='S';}
			else if(r.match[i]=='i'){r.match[i]='I';}
			else if(r.match[i]=='d'){r.match[i]='D';}
		}
	}

	/**
	 * Main implementation that fixes variants in a read by comparing against known variants.
	 * Processes the read's match string to identify substitutions and indels, then checks
	 * if these variants exist in the provided VarMap. Known variants have their match
	 * string characters modified to indicate their status.
	 *
	 * @param r Read to modify containing bases and match string
	 * @param sl SamLine containing alignment information for the read
	 * @param varMap Map of known variants for comparison
	 * @param scafMap Scaffold mapping information for coordinate resolution
	 * @return Number of variants successfully fixed in the read
	 */
	public static int fixVars(Read r, SamLine sl, VarMap varMap, ScafMap scafMap){
		if(r==null || r.bases==null || r.match==null){return 0;}
		assert(r.mapped());

		//		if(!Read.containsSubs(r.match)){return 0;}
		if(!Read.containsVars(r.match)){return 0;}
		final int scafnum=scafMap.getNumber(sl.rnameS());
		assert(scafnum>=0) : "Can't find scaffold "+sl.rnameS();
		if(scafnum<0){return 0;}
		//		System.err.println("A");

		if(r.match!=null && r.shortmatch()){
			r.toLongMatchString(false);
		}
		int varsFound=0;
		final byte[] match=r.match;
		final byte[] bases=r.bases;

		final boolean rcomp=(r.strand()==Shared.MINUS);
		if(rcomp){r.reverseComplementFast();}

		int rpos=sl.pos-1-SamLine.countLeadingClip(sl.cigar, true, true);

		//		System.err.println("varMap: \n"+varMap+"\n\n");
		byte prev='?';
		for(int bpos=0, mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
			assert(bpos<bases.length) : new String(match);
			final byte b=bases[bpos];

			if(m=='S'){
				Var v=new Var(scafnum, rpos, rpos+1, b, Var.SUB);
				if(varMap.containsKey(v)){
					varsFound++;
					match[mpos]='V';
					//					System.err.println("Found "+v+"\n");
					//				}else{
					//					System.err.println("Can't find "+v+" in\n"+varMap+"\n");
				}
			}else if(CallVariants.fixIndels && prev!=m && (m=='I' || m=='D')){
				int len=0;
				for(int i=mpos; i<match.length; i++){
					if(match[i]==m){len++;}else{break;}
				}
				byte replacement=Tools.toLowerCase(m);
				Var v;
				if(m=='D'){v=new Var(scafnum, rpos, rpos+len+1, 0, Var.DEL);}//Check the +1; may not be right
				else{
					byte[] alt=(len==1 ? Var.AL_MAP[b] : Arrays.copyOfRange(bases, bpos, bpos+len));
					v=new Var(scafnum, rpos, rpos, alt, Var.INS);
				}
				if(varMap.containsKey(v)){
					varsFound++;
					for(int i=mpos; i<match.length; i++){
						if(match[i]==m){match[i]=replacement;}else{break;}
					}
					//					System.err.println("Found "+v+"\n");
					//				}else{
					//					System.err.println("Can't find "+v+" in\n"+varMap+"\n");
				}
			}

			if(m!='D' && m!='d'){bpos++;}
			if(m!='I' && m!='i'){rpos++;}
			prev=m;
		}
		if(rcomp){r.reverseComplementFast();}

		//		assert(false) : new String(r.match);

		//		System.err.println("B:"+varsFound);
		return varsFound;
	}

	/**
	 * Finds unique substitution variants in a read that meet specified filtering criteria.
	 * Scans the read's alignment for substitutions and applies coverage, allele depth,
	 * allele fraction, and distance from read ends filters. Only substitutions are
	 * considered, not insertions or deletions.
	 *
	 * @param r Read to analyze for substitution variants
	 * @param sl SamLine containing alignment information
	 * @param varMap Map of known variants for comparison and filtering
	 * @param scafMap Scaffold mapping information for coordinate resolution
	 * @param maxVarDepth Maximum allowed variant depth for filtering
	 * @param maxAlleleFraction Maximum allowed allele fraction for filtering
	 * @param minCov Minimum coverage requirement for known variants
	 * @param minEDist Minimum distance from read ends to consider variants
	 * @return List of unique substitution variants meeting criteria, null if none found
	 */
	public static ArrayList<Var> findUniqueSubs(Read r, SamLine sl, VarMap varMap, ScafMap scafMap, int maxVarDepth, float maxAlleleFraction, int minCov, int minEDist){
		if(r==null || r.bases==null || r.match==null){return null;}
		assert(r.mapped());

		final int subs=Read.countSubs(r.match);
		if(subs==0){return null;}

		final int scafnum=scafMap.getNumber(sl.rnameS());
		assert(scafnum>=0) : "Can't find scaffold "+sl.rnameS();

		if(r.match!=null && r.shortmatch()){r.toLongMatchString(false);}

		final boolean rcomp=(r.strand()==Shared.MINUS);
		if(rcomp){r.reverseComplementFast();}

		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final ArrayList<Var> list=new ArrayList<Var>(subs);

		int rpos=sl.pos-1-SamLine.countLeadingClip(sl.cigar, true, true);
		int subsFound=0;
		for(int qpos=0, mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
			final byte b=bases[qpos];

			if(m=='S' && scafnum>=0){
				subsFound++;
				if(qpos>=minEDist && qpos<bases.length-minEDist){
					Var v=new Var(scafnum, rpos, rpos+1, b, Var.SUB);
					Var old=varMap.get(v);
					if(old==null){
						list.add(v);
					}else if(old.hasCoverage()){
						if(old.coverage()>=minCov){
							if(old.alleleCount()<=maxVarDepth || (maxAlleleFraction>0 && old.alleleFraction()<=maxAlleleFraction)){
								list.add(old);
							}
						}
					}else{
						if(old.alleleCount()<=maxVarDepth){list.add(old);}
					}
				}
			}

			if(m!='D'){qpos++;}
			if(m!='I'){rpos++;}
		}
		assert(subs==subsFound) : subs+", "+subsFound+", "+Read.countSubs(r.match)+"\n"+new String(match)+"\n"+new String(Read.toShortMatchString(r.match));
		if(rcomp){r.reverseComplementFast();}
		return list.isEmpty() ? null : list;
	}

	/**
	 * Finds unique variants of all types in a read that meet specified filtering criteria.
	 * More comprehensive than findUniqueSubs, handles substitutions, insertions, and
	 * deletions. Uses Var.toVars() to extract variants from the read alignment and
	 * applies filtering based on coverage and allele metrics.
	 *
	 * @param r Read to analyze for all variant types
	 * @param sl SamLine containing alignment information
	 * @param varMap Map of known variants for comparison and filtering
	 * @param scafMap Scaffold mapping information for coordinate resolution
	 * @param maxVarDepth Maximum allowed variant depth for filtering
	 * @param maxAlleleFraction Maximum allowed allele fraction for filtering
	 * @param minCov Minimum coverage requirement for known variants
	 * @param minEDist Minimum distance from read ends to consider variants
	 * @return List of unique variants meeting criteria, null if none found
	 */
	public static ArrayList<Var> findUniqueVars(Read r, SamLine sl, VarMap varMap, ScafMap scafMap, int maxVarDepth, float maxAlleleFraction, int minCov, int minEDist){
		if(r==null || r.bases==null || r.match==null){return null;}
		assert(r.mapped());

		final int vars=Read.countVars(r.match, Var.CALL_SUB, Var.CALL_INS, Var.CALL_DEL);
		if(vars==0){return null;}

		final int scafnum=scafMap.getNumber(sl.rnameS());
		assert(scafnum>=0) : "Can't find scaffold "+sl.rnameS();

		if(r.match!=null && r.shortmatch()){r.toLongMatchString(false);}

		final boolean rcomp=(r.strand()==Shared.MINUS);
		if(rcomp){
			r.reverseComplementFast();
			r.setSwapped(true);
		}

		final ArrayList<Var> list=Var.toVars(r, sl, false, scafMap);
		ArrayList<Var> list2=new ArrayList<Var>();
		for(Var v : list){
			if(v.endDistMax>=minEDist){
				Var old=varMap.get(v);
				if(old==null){
					list2.add(v);
				}else if(old.hasCoverage()){
					if(old.coverage()>=minCov){
						if(old.alleleCount()<=maxVarDepth || (maxAlleleFraction>0 && old.alleleFraction()<=maxAlleleFraction)){
							list2.add(old);
						}
					}
				}else{
					if(old.alleleCount()<=maxVarDepth){list2.add(old);}
				}
			}
		}
		if(rcomp){
			r.reverseComplementFast();
			r.setSwapped(false);
		}
		return list2.isEmpty() ? null : list2;
	}
	
	/**
	 * Loads variants from VCF files and marks them as forced variants.
	 * Processes comma-separated list of VCF file paths, loading all variants
	 * and adding them to the provided VarMap with forced status set to true.
	 * Forced variants are treated specially in variant calling pipelines.
	 *
	 * @param fnames Comma-separated list of VCF file paths to load
	 * @param scafMap Scaffold mapping information for variant positioning
	 * @param varMap Existing VarMap to add loaded variants to
	 * @param outstream PrintStream for progress and timing output
	 * @return The modified VarMap with added forced variants
	 */
	public static VarMap loadForcedVCF(String fnames, ScafMap scafMap, VarMap varMap, PrintStream outstream) {
		if(fnames==null){return null;}

		Timer t2=new Timer(outstream, true);
		String[] array=(fnames.indexOf(',')>=0 ? fnames.split(",") : new String[] {fnames});
		for(String fname : array){
			FileFormat ff=FileFormat.testInput(fname, FileFormat.VCF, null, true, false);
			VarMap varMap2=VcfLoader.loadFile(ff, scafMap, false, false);

			for(Var v : varMap2){
				v.clear();
				v.setForced(true);
				varMap.addUnsynchronized(v);
			}
		}

		t2.stop("Vars: \t"+varMap.size()+"\nTime: ");
		return varMap;
	}

}