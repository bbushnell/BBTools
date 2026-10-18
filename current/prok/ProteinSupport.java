package prok;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import dna.Data;
import dna.GeneticCode;
import prot.ProteinFamilyAssigner;
import prot.ProteinFamilyPacks;
import prot.HbmTailGate;
import prot.ProteinTailAssigner;
import prot.ProteinTailPool;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Optional immutable protein-family evidence for CallGenes, with thread-local search scratch.
 * Core selection never discovers additional installed packs. Functional labels do not
 * determine coding support; stop-facing sequence evidence does. Each generated stop
 * is searched once, before alternative starts are enumerated. Native PGM/NN inputs
 * remain unchanged; independent stop and ORF contributions enter afterward.
 * @author Keqing
 */
public final class ProteinSupport implements AutoCloseable {

	public static ProteinSupport load(final String mode, final String reference, final float bonus, final float stopBonus,
			final String[] outputs) throws IOException{
		return load(mode, reference, bonus, stopBonus, outputs, HbmTailGate.EXPERIMENTAL);
	}
	public static ProteinSupport load(final String mode, final String reference, final float bonus, final float stopBonus,
			final String[] outputs, final HbmTailGate.Policy policy) throws IOException{
		return load(mode, reference, bonus, stopBonus, outputs, policy, false);
	}
	public static ProteinSupport load(final String mode, final String reference, final float bonus, final float stopBonus,
			final String[] outputs, final HbmTailGate.Policy policy, final boolean modelWeighted) throws IOException{
		return load(mode, reference, bonus, stopBonus, outputs, policy, modelWeighted, Protection.AUTO);
	}
	public static ProteinSupport load(final String mode, final String reference, final float bonus, final float stopBonus,
			final String[] outputs, final HbmTailGate.Policy policy, final boolean modelWeighted, final Protection protection) throws IOException{
		return load(mode, reference, bonus, stopBonus, outputs, policy, modelWeighted, protection,
			searchPolicy("exhaustive", -1, 0, -1));
	}
	public static ProteinSupport load(final String mode, final String reference, final float bonus, final float stopBonus,
			final String[] outputs, final HbmTailGate.Policy policy, final boolean modelWeighted, final Protection protection,
			final ProteinTailAssigner.SearchPolicy searchPolicy) throws IOException{
		if(policy==null || protection==null || searchPolicy==null){throw new IllegalArgumentException("Protein support requires explicit tail, search and protection policies");}
		if(!("core".equals(mode) || "full".equals(mode)) || !Float.isFinite(bonus) || bonus<0
			|| !Float.isFinite(stopBonus) || stopBonus<0 || !Float.isFinite(bonus+stopBonus)){
			throw new IllegalArgumentException("Protein support requires annotation=core|full and finite nonnegative ORF/stop bonuses");
		}
		final String resolved=Data.findPath(reference, false);
		if(resolved==null || !Files.isDirectory(Paths.get(resolved))){
			throw new IllegalArgumentException("Missing BBGAP reference: "+reference+". Extract the BBGAP tar into BBTools/resources/.");
		}
		final Path root=Paths.get(resolved).toRealPath();
		final ProteinFamilyPacks.Pack[] packs=mode.equals("core") ? null : ProteinFamilyPacks.discover(root, false);
		final int n=packs==null ? 1 : packs.length;
		final ProteinFamilyAssigner[] libraries=new ProteinFamilyAssigner[n];
		final FamilyLabels[] labels=new FamilyLabels[n];
		final String[] ids=new String[n];
		for(int p=0; p<n; p++){
			final Path bundle=packs==null ? root : packs[p].root;
			final Path models=bundle.resolve("models").toRealPath();
			final Path compressed=bundle.resolve("labels/family_labels.tsv.gz");
			final Path label=Files.isRegularFile(compressed) ? compressed : bundle.resolve("labels/family_labels.tsv");
			protect(models, bundle, label, outputs);
			libraries[p]=packs==null ? ProteinFamilyAssigner.load(models) : packs[p].load();
			if(Files.isRegularFile(label)){labels[p]=FamilyLabels.load(label.toString(), libraries[p]);}
			ids[p]=packs==null ? "part1" : packs[p].id;
		}
		return new ProteinSupport(libraries, labels, ids, mode, bonus, stopBonus,
			packs==null ? libraries[0].modelPin() : ProteinFamilyPacks.referencePin(packs), policy, modelWeighted, protection, searchPolicy);
	}

	/** Check destinations before a caller writer can truncate a resource or hard-link alias. */
	private static void protect(final Path models, final Path bundle, final Path label,
			final String[] outputs) throws IOException{
		assert(outputs!=null) : "CallGenes must preflight all output streams before loading immutable protein resources";
		for(String name : outputs){
			if(!Tools.isOutputFileName(name)){continue;}
			final Path out=new File(name).getCanonicalFile().toPath();
			if(out.startsWith(bundle.toRealPath())){
				throw new IllegalArgumentException("Output would overwrite a BBGAP reference member: "+out);
			}
			if(!Files.exists(out)){continue;}
			try(DirectoryStream<Path> members=Files.newDirectoryStream(models)){
				for(Path member : members){rejectAlias(out, member);}
			}
			for(Path member : new Path[]{label, bundle.resolve("manifest.tsv"), bundle.resolve("manifest.tsv.gz"), bundle.resolve("pack_manifest.tsv")}){
				rejectAlias(out, member);
			}
		}
	}
	private static void rejectAlias(final Path output, final Path input) throws IOException{
		if(Files.exists(input) && Files.isSameFile(output, input)){
			throw new IllegalArgumentException("Output aliases a BBGAP reference member: "+output+" -> "+input);
		}
	}

	private ProteinSupport(ProteinFamilyAssigner[] libraries_, FamilyLabels[] labels_, String[] ids_,
			String mode_, float bonus_, float stopBonus_, String pin_, HbmTailGate.Policy policy_, boolean modelWeighted_, Protection protection_,
			ProteinTailAssigner.SearchPolicy searchPolicy_){
		assert(libraries_.length>0 && labels_.length==libraries_.length && ids_.length==libraries_.length)
			: "Every enabled model package needs a matching label slot and provenance ID";
		libraries=libraries_; labels=labels_; ids=ids_; mode=mode_; bonus=bonus_; stopBonus=stopBonus_; pin=pin_;
		modelWeighted=modelWeighted_;
		protection=protection_;
		long residues=0, families=0;
		for(ProteinFamilyAssigner library : libraries){
			for(int rank=0; rank<library.familyCount(); rank++){
				final int length=library.modelLength(rank);
				if(length<1){throw new IllegalArgumentException("Protein support requires positive model lengths");}
				residues=Math.addExact(residues, length); families=Math.addExact(families, 1);
			}
		}
		meanModelLength=residues/(double)families;
		middleBonuses=new float[libraries.length][];
		for(int p=0; p<libraries.length; p++){
			middleBonuses[p]=new float[libraries[p].familyCount()];
			for(int rank=0; rank<middleBonuses[p].length; rank++){
				final float addition=middleBonus(bonus, libraries[p].modelLength(rank), meanModelLength, modelWeighted);
				if(!Float.isFinite(addition+stopBonus)){throw new IllegalArgumentException("Model-scaled protein bonus overflows for family "+libraries[p].familyId(rank));}
				middleBonuses[p][rank]=addition;
			}
		}
		policy=policy_; searchPolicy=searchPolicy_; searchDensityText=Double.toString(searchPolicy.minSeedDensity);
		tails=new ProteinTailAssigner(policy, searchPolicy, libraries);
		pool=new ProteinTailPool(tails, Shared.threads());
	}
	public Worker newWorker(){return new Worker();}
	public int familyCount(){int count=0; for(ProteinFamilyAssigner library : libraries){count+=library.familyCount();} return count;}
	public String referencePin(){return pin;}
	/** Unweighted mean over every loaded family, never over observed hits or candidate lengths. */
	public double meanModelLength(){return meanModelLength;}
	/** -1 selects mode-specific defaults; exhaustive remains the frozen full-top50 comparator. */
	public static ProteinTailAssigner.SearchPolicy searchPolicy(String name, int minSeeds, double minSeedDensity, int lookahead){
		final ProteinTailAssigner.Mode mode;
		if("fast".equalsIgnoreCase(name)){mode=ProteinTailAssigner.Mode.FAST;}
		else if("exhaustive".equalsIgnoreCase(name)){mode=ProteinTailAssigner.Mode.EXHAUSTIVE;}
		else{throw new IllegalArgumentException("hbmsearch must be fast or exhaustive");}
		if(minSeeds<-1 || lookahead<-1){throw new IllegalArgumentException("HBM seed count and lookahead must be nonnegative or auto");}
		final boolean fast=mode==ProteinTailAssigner.Mode.FAST;
		return new ProteinTailAssigner.SearchPolicy(mode, minSeeds<0 ? (fast ? 1 : 0) : minSeeds,
			minSeedDensity, lookahead<0 ? (fast ? 4 : 0) : lookahead);
	}
	/** Controls only the final soft-pruning exemption; discovery and native start admission are untouched. */
	public enum Protection{
		AUTO("auto"), ON("on"), OFF("off");
		Protection(String name_){name=name_;}
		public static Protection forName(String name){
			for(Protection value : values()){if(value.name.equalsIgnoreCase(name)){return value;}}
			throw new IllegalArgumentException("hbmprotect must be auto, on, or off");
		}
		public boolean protects(float bonus, boolean supported){
			assert(Float.isFinite(bonus) && bonus>=0) : "ProteinSupport.addBonuses guarantees finite nonnegative evidence before pruning";
			return this==AUTO ? bonus>0 : this==ON && supported;
		}
		public final String name;
	}
	/** Annotation records exist for no-hit/unsupported ORFs too; only stop support can activate ON. */
	public static boolean protectsFromPruning(Orf orf){
		assert(orf!=null) : "Soft pruning operates on existing native candidates";
		final Evidence evidence=orf.proteinEvidence;
		// Preserve the old auto predicate for callers without attached protein metadata.
		return evidence==null ? orf.proteinBonus>0
			: evidence.reference.protection.protects(orf.proteinBonus, evidence.stop.supported());
	}
	/** Validated before any output is opened, including annotation=off invocations. */
	public static boolean modelWeighted(String name){
		if("flat".equalsIgnoreCase(name)){return false;}
		if("model".equalsIgnoreCase(name)){return true;}
		throw new IllegalArgumentException("hbmboostmode must be flat or model");
	}
	/** Scales only the middle term. No ORF/query length enters this shared-stop reward. */
	public static float middleBonus(float base, int length, double mean, boolean weighted){
		if(!Float.isFinite(base) || base<0 || length<1 || !Double.isFinite(mean) || mean<=0){
			throw new IllegalArgumentException("Middle bonus requires a finite nonnegative base and positive model/mean lengths");
		}
		final float value=weighted ? (float)(base*(length/mean)) : base;
		if(!Float.isFinite(value)){throw new IllegalArgumentException("Model-scaled middle bonus overflows");}
		return value;
	}
	@Override
	public void close(){pool.close();}

	/** ProteinTailAssigner searches top50 and may rerun each once with traceback. */
	public static final int DP_HISTOGRAM_BINS=101;
	/** One observation per searched query, never per alternative start or accepted family. */
	public static void recordDpCounts(final long[] matched, final long[] unmatched, final int fills, final boolean assigned){
		checkDpShape(matched); checkDpShape(unmatched);
		if(matched==unmatched || fills<0 || fills>=DP_HISTOGRAM_BINS){
			throw new IllegalArgumentException("DP query histograms must be distinct with fills in 0..100 (top50 plus trace reruns): "+fills);
		}
		final long[] counts=assigned ? matched : unmatched;
		if(counts[fills]<0){throw new IllegalArgumentException("Negative DP query count at fill bin "+fills);}
		counts[fills]=Math.addExact(counts[fills], 1);
	}
	/** Validate the entire merge before mutation so a bad donor cannot leave partial totals. */
	public static void mergeDpCounts(final long[] destination, final long[] source){
		checkDpShape(destination); checkDpShape(source);
		if(destination==source){throw new IllegalArgumentException("Cannot merge a DP histogram into itself");}
		for(int i=0; i<DP_HISTOGRAM_BINS; i++){
			if(destination[i]<0 || source[i]<0){throw new IllegalArgumentException("Negative DP query count at fill bin "+i);}
			Math.addExact(destination[i], source[i]);
		}
		for(int i=0; i<DP_HISTOGRAM_BINS; i++){destination[i]+=source[i];}
	}
	/** Literal bins include zero-work queries. Median averages the middle pair; p95
	 * uses nearest rank. Undefined empty statistics are NaN; empty maximum is zero.
	 * Aligned means at least one DP fill, not an accepted or positive-scoring match.
	 */
	public static String dpSummary(final long[] counts){
		checkDpShape(counts);
		long queries=0, fills=0;
		int max=0;
		for(int i=0; i<counts.length; i++){
			if(counts[i]<0){throw new IllegalArgumentException("Negative DP query count at fill bin "+i);}
			queries=Math.addExact(queries, counts[i]);
			fills=Math.addExact(fills, Math.multiplyExact(counts[i], (long)i));
			if(counts[i]>0){max=i;}
		}
		final long aligned=queries-counts[0];
		double median=Double.NaN, p95=Double.NaN;
		if(queries>0){
			final long lower=queries/2+queries%2, upper=queries/2+1, rank95=queries-queries/20;
			long cumulative=0;
			int low=-1, high=-1, percentile=-1;
			for(int i=0; i<counts.length; i++){
				cumulative+=counts[i]; // Total was checked above; every count is nonnegative.
				if(low<0 && cumulative>=lower){low=i;}
				if(high<0 && cumulative>=upper){high=i;}
				if(percentile<0 && cumulative>=rank95){percentile=i;}
			}
			assert(low>=0 && high>=0 && percentile>=0) : "Nonempty checked DP histogram must contain each requested rank";
			median=(low+high)*0.5; p95=percentile;
		}
		return "queries="+queries+" dp_fills="+fills+" zero_dp="+counts[0]
			+" mean="+(queries==0 ? Double.NaN : fills/(double)queries)+" median="+median+" p95="+p95+" max="+max
			+" aligned_queries="+aligned+" aligned_mean="+(aligned==0 ? Double.NaN : fills/(double)aligned);
	}
	private static void checkDpShape(final long[] counts){
		if(counts==null || counts.length!=DP_HISTOGRAM_BINS){
			throw new IllegalArgumentException("DP histogram requires exactly101 bins for top50 plus at most50 traceback reruns");
		}
	}
	/** Preserve legacy500AA GFF bytes; other scientific windows must be explicit. */
	public static void appendTailCap(final ByteBuilder row, final int cap){
		assert(row!=null && cap>0) : "Evidence retains a validated HbmTailGate policy and an existing GFF row";
		if(cap!=500){row.append(";bbgap_tail_cap=").append(cap);}
	}

	/** One instance per GeneCaller. Alignment failures must not leave stream consumers hanging. */
	public final class Worker {
		private Worker(){}
		/** Materialize immutable peptides before parallel work; no worker sees mutable genomic bases. */
		StopEvidence[] prepareStops(final ArrayList<Orf> longest, final byte[] orientedBases, final GeneticCode code){
			try{
				final GeneticCode effective=code==null ? GeneticCode.forTable(11) : code;
				final byte[][] queries=new byte[longest.size()][];
				final StopEvidence[] results=new StopEvidence[longest.size()];
				for(int i=0; i<longest.size(); i++){
					final Orf orf=longest.get(i);
					assert(orf.type==ProkObject.CDS && orf.scaflen==orientedBases.length)
						: "Stop-tail batches must preserve the contig coordinate frame";
					if(orf.stopTruncated){results[i]=new StopEvidence(null, "synthetic_stop", effective.id(), 0); continue;}
					final int first=orf.flipped()==0 ? orf.start : orientedBases.length-1-orf.stop;
					final int last=orf.flipped()==0 ? orf.stop : orientedBases.length-1-orf.start;
					final byte[] peptide=translateTail(orientedBases, first, last, effective, policy.maxLength);
					final String reason=translationReason(peptide);
					if(reason==null){queries[i]=peptide;}
					else{results[i]=new StopEvidence(null, reason, effective.id(), peptide.length);}
				}
				final ProteinTailAssigner.Hit[] hits=pool.assign(queries);
				// Joining publishes every fixed-index hit; only the caller updates counters.
				for(int i=0; i<results.length; i++){
					if(queries[i]!=null){
						results[i]=new StopEvidence(hits[i], null, effective.id(), queries[i].length);
						searched++; if(results[i].supported()){assigned++;}
						final ProteinTailAssigner.Hit hit=hits[i];
						recordDpCounts(matchedDpCounts, unmatchedDpCounts, hit.scored, hit.assigned());
						assert(!hit.compositionSkipped || (hit.shortlisted==0 && !hit.assigned() && hit.scored==0))
							: "Composition may be skipped only when the unchanged filters prove no candidate can be aligned";
						if(hit.compositionSkipped){compositionSkipped++;}
						shortlisted+=hit.shortlisted; lengthRejected+=hit.lengthRejected; seedRejected+=hit.seedRejected;
						visited+=hit.visited; scored+=hit.scored; recorded+=hit.recorded; accepted+=hit.accepted;
						unvisited+=hit.unvisited; if(hit.earlyStopped){earlyStopped++;}
					}else{unsupported++;}
				}
				return results;
			}catch(RuntimeException | Error failure){
				KillSwitch.kill("Protein stop batch failed before candidate scoring: "+failure);
				throw failure;
			}
		}
		void annotate(final Orf orf, final StopEvidence stop){
			try{attach(orf, stop);}catch(RuntimeException | AssertionError failure){
				KillSwitch.kill("Protein evidence attachment failed for "+orf.scafName+":"+orf.start+"-"+orf.stop+": "+failure);
				throw failure;
			}
		}
		private void attach(final Orf orf, final StopEvidence stop){
			assert(orf.type==ProkObject.CDS && orf.proteinEvidence==null && stop!=null)
				: "Each CDS candidate inherits its existing stop evidence exactly once";
			final float nativeScore=orf.orfScore;
			final float orfAddition=stop.supported() ? middleBonuses[stop.hit.packIndex][stop.hit.rank] : 0;
			final float stopAddition=stop.supported() ? stopBonus : 0;
			addBonuses(orf, stopAddition, orfAddition);
			orf.proteinEvidence=new Evidence(ProteinSupport.this, stop, stopAddition, orfAddition, nativeScore);
			if(orf.proteinBonus>0){boosted++; if(orf.proteinOnlyCandidate){rescuedStarts++;}}
		}
		public long searched, assigned, unsupported, withoutStarts, boosted, rescuedStarts;
		/** Counts accumulated after pool joins; scored counts actual DP matrix fills, including reruns. */
		public long shortlisted, lengthRejected, seedRejected, visited, scored, recorded, accepted, unvisited, earlyStopped;
		/** Queries whose search worker proved no candidate exists before computing composition. */
		public long compositionSkipped;
		/** Query counts by actual DP fills, partitioned by final assignment, including zero-work no-hits. */
		public final long[] matchedDpCounts=new long[DP_HISTOGRAM_BINS], unmatchedDpCounts=new long[DP_HISTOGRAM_BINS];
	}

	/** Complete stop is omitted. The arbitrary tail boundary never forces initiator methionine. */
	public static byte[] translateTail(byte[] bases, int start, int stop, GeneticCode code, int cap){
		if(bases==null || code==null || start<0 || stop>=bases.length || stop-start<2
			|| (stop-start+1)%3!=0 || cap<1){throw new IllegalArgumentException("Invalid biological stop-tail coordinates or cap");}
		if(!code.isStop(GeneticCode.codon(bases, stop-2))){throw new IllegalArgumentException("Tail endpoint is not a stop under code "+code.id());}
		final int count=Math.min(cap, (stop-start+1)/3-1), first=stop-2-3*count;
		final byte[] peptide=new byte[count];
		for(int i=0; i<count; i++){peptide[i]=code.aminoAcid(GeneticCode.codon(bases, first+3*i));}
		return peptide;
	}
	/** Adds evidence once after native NN/minor-start scoring, without altering those inputs. */
	public static void addBonuses(Orf orf, float stopAddition, float orfAddition){
		final float total=stopAddition+orfAddition;
		if(orf==null || orf.proteinBonus!=0 || !Float.isFinite(total) || stopAddition<0 || orfAddition<0){
			throw new IllegalArgumentException("Protein score additions must be finite, nonnegative and applied once");
		}
		orf.proteinBonus=total;
		if(total>0){orf.orfScore+=total;}
	}
	/** Shared by every alternative start at this stop; no candidate-dependent search result. */
	public static final class StopEvidence{
		public StopEvidence(ProteinTailAssigner.Hit hit_, String reason_, int table_, int length_){
			hit=hit_; reason=reason_; table=table_; length=length_;
		}
		public boolean supported(){return hit!=null && hit.assigned();}
		final ProteinTailAssigner.Hit hit;
		final String reason;
		final int table, length;
	}

	public static String translationReason(final byte[] peptide){
		if(peptide==null || peptide.length==0){return "empty_translation";}
		for(byte residue : peptide){
			if(residue=='*'){return "internal_stop";}
			if(residue=='X'){return "ambiguous_codon";}
		}
		return null;
	}

	/** Detached evidence survives candidate sorting and DP; score terms remain separable. */
	static final class Evidence {
		Evidence(ProteinSupport reference_, StopEvidence stop_, float stopAddition_, float orfAddition_, float nativeScore_){
			assert(reference_!=null && stop_!=null && Float.isFinite(stopAddition_+orfAddition_))
				: "Every attached score bonus must retain its immutable reference and a finite nonnegative magnitude";
			reference=reference_; stop=stop_; stopAddition=stopAddition_; orfAddition=orfAddition_; nativeScore=nativeScore_;
		}
		void append(final ByteBuilder row){
			assert(row!=null) : "Protein annotations append to, never replace, existing CallGenes attributes";
			final ProteinTailAssigner.Hit hit=stop.hit;
			final boolean assigned=stop.supported();
			final int p=assigned ? hit.packIndex : 0;
			row.append(";bbgap_status=").append(stop.reason!=null ? "unsupported" : assigned ? "stop_supported" : "no_hit");
			row.append(";bbgap_mode=").append(reference.mode).append(";bbgap_reference=").append(reference.pin);
			row.append(";bbgap_rule=").append(reference.tails.rule()).append(";bbgap_translation_table=").append(stop.table);
			final boolean fast=reference.searchPolicy.mode==ProteinTailAssigner.Mode.FAST;
			if(fast){
				row.append(";bbgap_search=fast;bbgap_competitor_scope=visited;bbgap_min_seeds=").append(reference.searchPolicy.minSeeds);
				row.append(";bbgap_min_seed_density=").append(reference.searchDensityText);
				row.append(";bbgap_lookahead=").append(reference.searchPolicy.lookahead);
			}
			row.append(";bbgap_evidence_scope=stop_tail;bbgap_tail_length=").append(stop.length);
			appendTailCap(row, reference.policy.maxLength);
			row.append(";bbgap_native_score=").append(nativeScore, 6).append(";bbgap_call_bonus=").append(stopAddition+orfAddition, 6);
			row.append(";bbgap_stop_bonus=").append(stopAddition, 6).append(";bbgap_orf_bonus=").append(orfAddition, 6);
			if(reference.protection!=Protection.AUTO){
				row.append(";bbgap_protection=").append(reference.protection.name);
				row.append(";bbgap_soft_pruning_exempt=").append(reference.protection.protects(stopAddition+orfAddition, assigned) ? 1 : 0);
			}
			if(reference.modelWeighted){
				row.append(";bbgap_bonus_rule=model_length_v1;bbgap_mean_model_length=").append(reference.meanModelLength, 8);
				if(assigned){row.append(";bbgap_model_length=").append(reference.libraries[p].modelLength(hit.rank));}
			}
			if(stop.reason!=null){row.append(";bbgap_reason=").append(stop.reason);}
			if(!assigned){return;}
			row.append(";bbgap_pack=").append(reference.ids[p]).append(";bbgap_package=").append(reference.libraries[p].modelPin());
			row.append(";bbgap_family_id=").append(hit.familyId).append(";bbgap_rep_id="); AnnotationGff.escaped(row, hit.representative);
			row.append(";bbgap_score64=").append(hit.score64).append(";bbgap_identity=").append(hit.identity, 8);
			row.append(";bbgap_paired_fraction=").append(hit.pairedFraction, 8).append(";bbgap_core_cov=").append(hit.coreCoverage, 8);
			row.append(";bbgap_tail_paired=").append(hit.paired).append(";bbgap_tail_terminal_paired=").append(hit.terminalPaired);
			row.append(";bbgap_tail_query_suffix=").append(hit.querySuffix).append(";bbgap_tail_model_suffix=").append(hit.referenceSuffix);
			row.append(";bbgap_tail_query_start=").append(hit.queryStart).append(";bbgap_tail_model_start=").append(hit.referenceStart);
			if(reference.policy.profileEvidence){
				row.append(";bbgap_identity_scope=paired;bbgap_tail_terminal_score64=").append(hit.terminalScore64);
			}
			if(fast){
				row.append(";bbgap_visited_runner_score64=").append(hit.runnerUpScore64).append(";bbgap_visited_other_score64=").append(hit.bestOtherScore64);
			}else{row.append(";bbgap_runner_supported_score64=").append(hit.runnerUpScore64).append(";bbgap_best_other_score64=").append(hit.bestOtherScore64);}
			final FamilyLabels label=reference.labels[p];
			if(label!=null){row.append(";bbgap_label_version="); AnnotationGff.escaped(row, label.version);}
			AnnotateGenes.appendFamilyLabel(row, label, hit.rank);
		}
		final ProteinSupport reference;
		final StopEvidence stop;
		final float stopAddition, orfAddition, nativeScore;
	}
	private final ProteinFamilyAssigner[] libraries;
	private final FamilyLabels[] labels;
	private final String[] ids;
	private final ProteinTailAssigner tails;
	private final HbmTailGate.Policy policy;
	private final ProteinTailAssigner.SearchPolicy searchPolicy;
	private final Protection protection;
	private final ProteinTailPool pool;
	private final String mode, pin, searchDensityText;
	private final float bonus, stopBonus;
	private final boolean modelWeighted;
	private final double meanModelLength;
	private final float[][] middleBonuses;
}
