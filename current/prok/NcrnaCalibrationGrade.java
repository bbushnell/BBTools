package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import map.ObjectSet;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

import prok.NcrnaGffFamilyGrader.FamilyMap;
import prok.NcrnaGffFamilyGrader.GenericGradeResult;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;

/** One-to-one historical/full-core grades with explicit unresolved-call accounting.
 * Exports apparent FP loci for independent adjudication; never declares them false.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCalibrationGrade {
	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, NcrnaCalibrationGrade.class, false);
		String calls=null, source=null, full=null, unresolved=null, map=null, out=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){calls=value;}else if(key.equals("source")){source=value;}else if(key.equals("full")){full=value;}
			else if(key.equals("unresolved")){unresolved=value;}else if(key.equals("familymap")){map=value;}
			else if(key.equals("outprefix")){out=value;}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(calls!=null && source!=null && full!=null && unresolved!=null && map!=null && out!=null,
			"Bind actual caller GFF, both coordinate views, unresolved masks and family attribution");
		final FamilyMap families=NcrnaGffFamilyGrader.loadFamilyMap(map);
		require(families.prefixByFamily.size()==1, "This development report grades one explicitly selected family at a time");
		final ArrayList<Locus> query=NcrnaGffFamilyGrader.loadCallsFromGff(calls, families);
		final ArrayList<Locus> original=NcrnaGffFamilyGrader.loadTruthTsv(source, families), cores=NcrnaGffFamilyGrader.loadTruthTsv(full, families);
		final ArrayList<Locus> masks=NcrnaGffFamilyGrader.loadTruthTsv(unresolved, families);
		require(original.size()==cores.size()+masks.size(), "Supported and unresolved coordinate views must conserve the original locus population");
		RrnaResourceIO.fresh(out, ".source.summary.tsv", ".source.matches.tsv", ".full.summary.tsv", ".full.matches.tsv",
			".calls.tsv", ".novel.gff", ".novel.contigs.txt", ".accounting.tsv");
		final ArrayList<Locus> unique=new ArrayList<Locus>(), scored=new ArrayList<Locus>();
		final IntList uniqueIndices=new IntList(), indices=new IntList();final boolean[] unscored=new boolean[query.size()], duplicate=new boolean[query.size()];
		final ObjectSet<String> physical=new ObjectSet<String>(String.class);
		int duplicates=0;
		for(int i=0; i<query.size(); i++){
			final Locus c=query.get(i);boolean fullEdge=false, unknownOverlap=false;
			duplicate[i]=!physical.add(c.seqid+'\t'+c.family+'\t'+c.strand+'\t'+c.start+'\t'+c.stop);
			if(duplicate[i]){duplicates++;continue;}
			unique.add(c);uniqueIndices.add(i);
			for(Locus t:cores){if(NcrnaGffFamilyGrader.overlapsReciprocal(c, t, .5f)){fullEdge=true;break;}}
			if(!fullEdge){for(Locus t:masks){if(c.seqid.equals(t.seqid) && c.strand==t.strand && c.start<=t.stop && c.stop>=t.start){unknownOverlap=true;break;}}}
			unscored[i]=unknownOverlap;
			if(!unknownOverlap){scored.add(c);indices.add(i);}
		}
		final MatchResult sourceMatch=NcrnaGffFamilyGrader.matchOneToOne(unique, original, .5f), coreMatch=NcrnaGffFamilyGrader.matchOneToOne(scored, cores, .5f);
		final boolean[] sourceTp=new boolean[query.size()], coreTp=new boolean[query.size()];
		for(int i=0; i<uniqueIndices.size; i++){sourceTp[uniqueIndices.get(i)]=sourceMatch.callConsumed[i];}
		for(int i=0; i<indices.size; i++){coreTp[indices.get(i)]=coreMatch.callConsumed[i];}
		final GenericGradeResult sourceGrade=NcrnaGffFamilyGrader.gradeGeneric(original, unique, Collections.emptySet(), families);
		final GenericGradeResult fullGrade=NcrnaGffFamilyGrader.gradeGeneric(cores, scored, Collections.emptySet(), families);
		for(int i=0; i<query.size(); i++){
			if(!duplicate[i]){continue;}
			final String family=query.get(i).family;
			sourceGrade.counts.get(family).fp++;sourceGrade.counts.get(family).dup++;sourceGrade.unionFP++;
			fullGrade.counts.get(family).fp++;fullGrade.counts.get(family).dup++;fullGrade.unionFP++;
		}
		NcrnaGffFamilyGrader.writeGenericReport(out+".source", sourceGrade);
		NcrnaGffFamilyGrader.writeGenericReport(out+".full", fullGrade);
		write(out, query, sourceTp, coreTp, unscored, duplicate, scored.size()+duplicates);Shared.closeStream(pp.outstream);
	}
	private static void write(final String out, final ArrayList<Locus> query, final boolean[] sourceTp, final boolean[] coreTp,
		final boolean[] unscored, final boolean[] duplicate, final int scored){
		assert(query.size()==coreTp.length && query.size()==unscored.length):"Every emitted call needs a disposition in both coordinate views";
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".calls.tsv"), novel=RrnaResourceIO.writer(out+".novel.gff"), contigs=RrnaResourceIO.writer(out+".novel.contigs.txt");
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ArrayList<String> names=new ArrayList<String>();
		final ByteBuilder b=new ByteBuilder();int unknown=0, candidates=0;
		try{
			rows.println("call_id\tseqid\tfamily\tstart1\tstop1\tstrand\tsource_outcome\tfull_outcome\texact_duplicate\tindependent_review");novel.println("##gff-version 3");
			for(int i=0; i<query.size(); i++){
				final Locus c=query.get(i);if(unscored[i]){unknown++;}
				final String full=unscored[i]?"UNSCORED_UNRESOLVED_OVERLAP":coreTp[i]?"TP":"APPARENT_FP";
				final boolean review=!duplicate[i] && (!sourceTp[i] || !unscored[i] && !coreTp[i]);
				rows.print(b.clear().append("candidate_").append(i).tab().append(c.seqid).tab().append(c.family).tab().append(c.start).tab().append(c.stop)
					.tab().append(c.strand).tab().append(sourceTp[i]?"TP":"APPARENT_FP").tab().append(full).tab().append(duplicate[i]).tab().append(review).nl());
				if(review){
					candidates++;
					novel.print(b.clear().append(c.seqid).append("\tCalibration\tRNA\t").append(c.start).tab().append(c.stop)
						.append("\t.\t").append(c.strand).append("\t.\tID=candidate_").append(i).nl());
					if(seen.add(c.seqid)){names.add(c.seqid);}
				}
			}
			Collections.sort(names);for(String name:names){contigs.println(name);}
		}finally{RrnaResourceIO.finish(rows, novel, contigs);}
		require(scored+unknown==query.size(), "Full-view scored/unscored calls must conserve all emitted calls");
		final ByteStreamWriter accounting=RrnaResourceIO.writer(out+".accounting.tsv");
		try{
			accounting.println("calls\tfull_scored_calls\tunscored_calls\tindependent_review_candidates");
			accounting.print(b.clear().append(query.size()).tab().append(scored).tab().append(unknown).tab().append(candidates).nl());
		}finally{RrnaResourceIO.finish(accounting);}
		System.err.println("NCRNA_CALIBRATION_GRADE_PASS calls="+query.size()+" unscored="+unknown+" independent_review="+candidates+"; apparent counts, not adjudicated precision");
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
}
