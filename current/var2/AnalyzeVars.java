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
 * Marks known variants in read match strings and finds variants absent or
 * weakly represented in an existing VarMap.
 * <p>
 * This is used by BBDuk and CalcTrueQuality to distinguish already-known
 * substitutions, insertions, and deletions from variants that remain useful
 * for masking or calibration. VCF loading here is limited to loading forced
 * variants into an existing VarMap.
 *
 * @author Brian Bushnell
 * @author Isla
 * @date December 2024
 */
public class AnalyzeVars {

/*--------------------------------------------------------------*/
/*----------------      Known-Variant Marking   ----------------*/
/*--------------------------------------------------------------*/

	/**
	 * Marks variants in a mapped read's match string when they are already in
	 * varMap.
	 * <p>
	 * This wrapper rejects null, unmapped, and unaligned reads before delegating
	 * to the alignment-aware implementation.
	 *
	 * @param r Read to modify
	 * @param varMap Known variants
	 * @param scafMap Scaffold names to numeric IDs
	 * @return Number of known variants marked
	 */
	public static int fixVars(Read r, VarMap varMap, ScafMap scafMap){
		if(r==null || r.bases==null || r.match==null || r.samline==null){return 0;}
		final SamLine sl=r.samline;
		if(!sl.mapped()){return 0;}
		return AnalyzeVars.fixVars(r, sl, varMap, scafMap);
	}

	/**
	 * Restores match symbols changed by {@link #fixVars(Read, VarMap, ScafMap)}.
	 *
	 * @param r Read to restore; null-safe
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
	 * Marks known substitutions as 'V', known insertions as 'i', and known
	 * deletions as 'd' in a mapped read's long match string.
	 * <p>
	 * Minus-strand reads are reverse-complemented temporarily so match-string
	 * coordinates can be compared directly to VCF-derived VarMap coordinates.
	 *
	 * @param r Read to modify
	 * @param sl Alignment record for the read
	 * @param varMap Known variants
	 * @param scafMap Scaffold names to numeric IDs
	 * @return Number of known variants marked
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
				//Checked (G11 2026-06-29): the +1 is CORRECT. The varMap here is loaded via VcfLoader->
				//VcfToVar.fromVCF, which uses anchor-INCLUSIVE DEL coords (stop=pos+len+1); this rpos+len+1
				//matches it, so known deletions DO get found. (It does NOT match the anchor-exclusive
				//toSubsAndIndels/STA-STO convention, but that is not the parser this compares against.)
				if(m=='D'){v=new Var(scafnum, rpos, rpos+len+1, 0, Var.DEL);}
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

/*--------------------------------------------------------------*/
/*----------------       Unique-Variant Calls   ----------------*/
/*--------------------------------------------------------------*/

	/**
	 * Finds substitution calls that are absent from varMap or weak enough in
	 * varMap to keep for masking.
	 *
	 * @param r Read to analyze for substitution variants
	 * @param sl Alignment record for the read
	 * @param varMap Known variants and their observed depths
	 * @param scafMap Scaffold names to numeric IDs
	 * @param maxVarDepth Maximum known-variant allele depth to keep
	 * @param maxAlleleFraction Maximum known-variant allele fraction to keep
	 * @param minCov Minimum coverage before allele-depth filters apply
	 * @param minEDist Minimum end distance for emitted variants
	 * @return Unique substitution calls, or null if none pass
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
	 * Finds substitution, insertion, and deletion calls that are absent from
	 * varMap or weak enough in varMap to keep for masking.
	 *
	 * @param r Read to analyze for all variant types
	 * @param sl Alignment record for the read
	 * @param varMap Known variants and their observed depths
	 * @param scafMap Scaffold names to numeric IDs
	 * @param maxVarDepth Maximum known-variant allele depth to keep
	 * @param maxAlleleFraction Maximum known-variant allele fraction to keep
	 * @param minCov Minimum coverage before allele-depth filters apply
	 * @param minEDist Minimum end distance for emitted variants
	 * @return Unique variant calls, or null if none pass
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

/*--------------------------------------------------------------*/
/*----------------            Loading           ----------------*/
/*--------------------------------------------------------------*/
	
	/**
	 * Loads comma-delimited VCF paths, clears each loaded variant's counters,
	 * marks it forced, and adds it to the supplied map.
	 *
	 * @param fnames Comma-delimited VCF paths to load
	 * @param scafMap Scaffold names to numeric IDs
	 * @param varMap Map to receive forced variants
	 * @param outstream Progress stream
	 * @return The modified map, or null if {@code fnames} is null
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
