package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;

/** Saved-table profiling for BuildNcrnaReference, with no genomic payload reads.
 * Span quantiles weight raw CM-hit events, not genes, and do not choose cutoffs.
 * Availability means a GFF was parsed; even a parsed empty GFF is not a negative
 * biological reference. Every selected input's original path claims are conserved.
 * Instances and all accumulation are confined to the calling thread.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceProfiles{

	static void process(final String input, final String output, final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceProfiles profiles=load(input, catalog, roster, true);
		profiles.writeModels(output+".model_profile.tsv");
		profiles.writeAvailability(output);
	}

	/** Loads saved population/count bindings without writing another profile. */
	static NcrnaReferenceProfiles load(final String input, final NcrnaEvidenceCatalog catalog, final Roster roster){
		return load(input, catalog, roster, false);
	}
	/** Batch-aware consumers opt in explicitly; unadapted decision tools fail before use. */
	static NcrnaReferenceProfiles load(final String input, final NcrnaEvidenceCatalog catalog, final Roster roster, final boolean allowBatch){
		assert(catalog!=null && roster!=null) : "Profiles join saved rows to the original census and queried-model roster";
		final NcrnaReferenceProfiles profiles=new NcrnaReferenceProfiles(catalog, roster);
		profiles.readSummary(input+".summary.tsv");
		require(allowBatch || !profiles.batchEvidence, "This consumer does not support the captured batch evidence schema");
		profiles.readCoverage(input+".coverage.tsv");
		profiles.readSources(input+".sources.tsv.gz");
		return profiles;
	}
	/** Source-capture binding does not depend on which CM roster was searched.
	 * This view validates summary/source claims only; it has no model histogram. */
	static NcrnaReferenceProfiles sourceClaims(final String input, final NcrnaEvidenceCatalog catalog){
		final NcrnaReferenceProfiles profiles=new NcrnaReferenceProfiles(catalog, null);
		profiles.readSummary(input+".summary.tsv");
		profiles.readSources(input+".sources.tsv.gz");
		return profiles;
	}
	int observed(){return observed;}
	long events(){return hits;}
	long modelEvents(final int index){return models.get(index).sum;}
	long modelPairs(){return modelPairs;}
	long historicalPairs(){return historicalPairs;}
	String evidenceState(){return evidenceState;}
	boolean batchEvidence(){return batchEvidence;}
	boolean selected(final String id){return byInput.get(id).seenClaims>0;}
	long retained(){return retained;}
	long retained(final String id){return byInput.get(id).retained;}
	/** Compares canonical saved source claims by their catalog slots, independently
	 * of worker/output order. These small metadata rows retain exact paths, states,
	 * membership and feature accounting; they do not certify live GFF bytes. */
	boolean sameSources(final NcrnaReferenceProfiles other, final String id){
		final Availability a=byInput.get(id), b=other.byInput.get(id);
		if(a==null || b==null || a.claims!=b.claims || a.seenClaims!=b.seenClaims){return false;}
		for(int i=0; i<a.claims; i++){
			if(!Arrays.equals(a.sourceRows[i], b.sourceRows[i])){return false;}
		}
		return true;
	}
	/** Parsed claims, including empty GFFs; availability never establishes gene absence. */
	int availableClaims(final String id){
		final Availability input=byInput.get(id);
		return input.nonempty+input.empty;
	}

	private NcrnaReferenceProfiles(final NcrnaEvidenceCatalog catalog_, final Roster roster_){
		catalog=catalog_;
		roster=roster_;
		for(int i=0; roster!=null && i<roster.queried(); i++){
			final ModelProfile model=new ModelProfile(roster.query(i));
			models.add(model);
			byModel.put(model.model.accession, model);
		}
		final ByteBuilder buffer=new ByteBuilder();
		for(int i=0; i<catalog.size(); i++){
			final Availability item=new Availability(catalog.get(i), buffer);
			inputs.add(item);
			byInput.put(item.genome.id, item);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved inputs         ----------------*/
	/*--------------------------------------------------------------*/

	private void readSummary(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.SUMMARY_HEADER, path);
			final byte[] row=in.reader.nextLine();
			require(row!=null, "Missing saved integration summary");
			final LineParser1 p=new LineParser1('\t');
			p.set(row);
			require(p.terms()==16, "Saved integration summary needs sixteen fields");
			observed=p.parseInt(0);
			bases=p.parseLong(2);
			contigs=p.parseLong(3);
			hits=p.parseLong(4);
			claims=p.parseLong(5);
			unavailable=p.parseLong(6);
			features=p.parseLong(7);
			retained=p.parseLong(8);
			excluded=p.parseLong(9);
			historicalPairs=p.parseLong(10);
			modelPairs=p.parseLong(11);
			final int queried=p.parseInt(12);
			require(observed>0 && observed<=catalog.size() && p.parseInt(1)==catalog.size()
				&& queried>0 && (roster==null || queried==roster.queried()), "Saved population/model denominators disagree with supplied metadata");
			require(hits>=0 && claims>=observed && unavailable>=0 && unavailable<=claims
				&& retained>=0 && excluded>=0 && features==Math.addExact(retained, excluded), "Invalid saved source/event accounting");
			require(p.parseLong(13)==(long)observed*queried, "Every observed input needs all queried-model completions");
			require(bases>0 && contigs>0 && p.parseLong(10)>=0 && p.parseLong(11)>=0
				&& p.parseLong(14)>=0 && p.parseLong(14)<=p.parseLong(13), "Invalid saved dimensions, overlap counts or zero-hit pairs");
			batchEvidence=p.termEquals(BuildNcrnaReference.BATCH_STATUS, 15);
			require(!batchEvidence || observed==catalog.size(), "Batch evidence must preserve its complete supplied population");
			evidenceState=batchEvidence ? BuildNcrnaReference.BATCH_STATUS : observed==catalog.size() ? "EVIDENCE_COMPLETE_NOT_ADJUDICATED" : "SAMPLE_EVIDENCE_NOT_ADJUDICATED";
			require(p.termEquals(evidenceState, 15) && in.reader.nextLine()==null, "Expected one completed evidence summary with the proper population scope");
		}
	}

	private void readCoverage(final String path){
		final LineParser1 p=new LineParser1('\t');
		ModelProfile current=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.COVERAGE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==8, "Saved coverage rows need eight fields");
				if(current==null || !p.termEquals(current.model.accession, 0)){current=byModel.get(text(p, 0));}
				require(current!=null, "Coverage uses an unqueried or unknown accession");
				require(p.termEquals(current.model.release, 1) && p.parseInt(2)==current.model.length && p.parseInt(7)==observed,
					"Coverage release, model length and input count must match the supplied roster/summary");
				final long count=p.parseLong(5), denominator=p.parseLong(6);
				require(denominator>=0 && (current.denominator<0 || current.denominator==denominator), "Inconsistent accession hit denominator");
				current.denominator=denominator;
				if(p.termEquals("NA", 3)){
					require(p.termEquals("NA", 4) && count==0 && denominator==0 && !current.zeroRow,
						"Only one explicit zero-hit row may use NA span/truncation");
					current.zeroRow=true;
				}else{
					final int span=p.parseInt(3), truncation=p.parseInt(4);
					require(span>0 && span<=current.model.length && truncation>=0 && truncation<4 && count>0 && denominator>0,
						"Numeric coverage bins require valid spans/truncation and positive event counts");
					final int cell=span*4+truncation;
					require(current.bins[cell]==0, "Repeated accession/span/truncation bin");
					current.bins[cell]=count;
					current.sum=Math.addExact(current.sum, count);
					current.truncations[truncation]=Math.addExact(current.truncations[truncation], count);
				}
			}
		}
		long total=0;
		for(final ModelProfile model:models){
			require(model.denominator>=0 && model.sum==model.denominator, "Missing queried-model profile or unconserved hit histogram");
			total=Math.addExact(total, model.sum);
		}
		require(total==hits, "All queried-model histograms must conserve the saved raw-event total");
	}

	private void readSources(final String path){
		final LineParser1 p=new LineParser1('\t');
		Availability current=null;
		long totalClaims=0, totalUnavailable=0, totalFeatures=0, totalRetained=0, totalExcluded=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.SOURCE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==12, "Saved source rows need twelve fields");
				if(current==null || !p.termEquals(current.genome.id, 0)){current=byInput.get(text(p, 0));}
				require(current!=null, "Saved annotation claim names an unknown catalog input");
				final int claim=current.find(p);
				require(claim>=0 && !current.seen[claim] && p.parseInt(3)==current.members[claim],
					"Source claim is unknown, repeated or has a different membership count");
				current.seen[claim]=true;
				current.sourceRows[claim]=row;
				current.seenClaims++;
				final boolean parsed=p.termEquals("PARSED_GFF", 2), empty=p.termEquals("EMPTY_GFF", 2);
				final long sourceRows=p.parseLong(4), kept=p.parseLong(5), dropped=p.parseLong(6);
				require(sourceRows>=0 && kept>=0 && dropped>=0 && sourceRows==Math.addExact(kept, dropped), "Unconserved historical source rows");
				long exclusions=0;
				for(int field=7; field<=9; field++){
					final long count=p.parseLong(field);
					require(count>=0, "Negative historical exclusion count");
					exclusions=Math.addExact(exclusions, count);
				}
				require(exclusions==dropped && p.parseLong(10)>=0 && p.parseLong(11)>=0, "Invalid exclusion/structural accounting");
				if(parsed){
					require(sourceRows>0, "A parsed nonempty GFF must have source features");
					current.nonempty++;
				}else if(empty){
					require(sourceRows==0, "An empty GFF cannot have source features");
					current.empty++;
				}else{
					boolean known=false;
					for(final String state:UNAVAILABLE){known|=p.termEquals(state, 2);}
					require(known && sourceRows==0 && p.parseLong(10)==0 && p.parseLong(11)==0,
						"Unavailable GFF state must be explicit and must not invent parsed features");
					current.unavailable++;
					totalUnavailable++;
				}
				current.features=Math.addExact(current.features, sourceRows);
				current.retained=Math.addExact(current.retained, kept);
				current.excluded=Math.addExact(current.excluded, dropped);
				totalClaims++;
				totalFeatures=Math.addExact(totalFeatures, sourceRows);
				totalRetained=Math.addExact(totalRetained, kept);
				totalExcluded=Math.addExact(totalExcluded, dropped);
			}
		}
		int selected=0;
		long selectedBases=0, selectedContigs=0;
		for(final Availability item:inputs){
			if(item.seenClaims==0){continue;}
			require(item.seenClaims==item.claims, "Every selected input must retain all original GFF claims");
			selected++;
			selectedBases=Math.addExact(selectedBases, item.genome.bases);
			selectedContigs=Math.addExact(selectedContigs, item.genome.records);
		}
		require(selectedBases==bases && selectedContigs==contigs, "Saved input dimensions disagree with the selected catalog population");
		require(selected==observed && totalClaims==claims && totalUnavailable==unavailable && totalFeatures==features
			&& totalRetained==retained && totalExcluded==excluded, "Source profiles disagree with saved population/accounting totals");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reports            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeModels(final String path){
		try(Output out=new Output(path)){
			out.writer.println("accession\tbiological_label\treference_family\tevidence_role\trelease\tmodel_length\tinputs_observed\traw_hit_events\tuntruncated_events\tnative_5prime_truncated\tnative_3prime_truncated\tnative_both_truncated\tuntruncated_full_model_span\ttruncated_full_model_span\tp05_model_span\tmedian_model_span\tp95_model_span\tmodal_model_span\tscoring_scope");
			final ByteBuilder b=new ByteBuilder();
			for(final ModelProfile profile:models){
				final Model model=profile.model;
				BuildNcrnaReference.field(b.clear(), model.accession).tab();
				BuildNcrnaReference.field(b, model.scope.label).tab();
				BuildNcrnaReference.field(b, model.scope.family).tab();
				BuildNcrnaReference.field(b, model.scope.role).tab();
				BuildNcrnaReference.field(b, model.release).tab().append(model.length).tab().append(observed).tab().append(profile.sum).tab();
				for(final long count:profile.truncations){b.append(count).tab();}
				final int end=model.length*4;
				b.append(profile.bins[end]).tab().append(profile.bins[end+1]+profile.bins[end+2]+profile.bins[end+3]).tab();
				if(profile.sum==0){b.append("NA\tNA\tNA\tNA\t");}
				else{
					b.append(profile.quantile(5)).tab().append(profile.quantile(50)).tab().append(profile.quantile(95)).tab().append(profile.mode()).tab();
				}
				BuildNcrnaReference.field(b, model.scope.scoringScope).nl();
				out.writer.print(b);
			}
		}
	}

	private void writeAvailability(final String prefix){
		final long[][] totals=new long[3][6];
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+".input_availability.tsv")){
			out.writer.println("input_id\tavailable_claims\tunavailable_claims\tnonempty_parsed_claims\tempty_parsed_claims\tsource_feature_rows\tretained_rows\texcluded_rows\tstatus");
			for(final Availability item:inputs){
				if(item.seenClaims==0){continue;}
				final int available=item.nonempty+item.empty;
				final int category=available==0 ? 0 : item.unavailable==0 ? 1 : 2;
				final long[] row=totals[category];
				row[0]++;
				row[1]+=item.claims;
				row[2]+=available;
				row[3]+=item.unavailable;
				if(available>0 && item.nonempty==0){row[4]++;}
				if(item.nonempty>0){row[5]++;}
				BuildNcrnaReference.field(b.clear(), item.genome.id).tab().append(available).tab().append(item.unavailable).tab();
				b.append(item.nonempty).tab().append(item.empty).tab().append(item.features).tab().append(item.retained).tab().append(item.excluded).tab();
				out.writer.print(b.append(CATEGORIES[category]).nl());
			}
		}
		try(Output out=new Output(prefix+".availability.tsv")){
			out.writer.println("category\tinputs\tsource_claims\tavailable_claims\tunavailable_claims\tinputs_with_only_empty_available_gff\tinputs_with_nonempty_parsed_gff");
			for(int i=0; i<totals.length; i++){
				b.clear().append(CATEGORIES[i]);
				for(final long count:totals[i]){b.tab().append(count);}
				out.writer.print(b.nl());
			}
		}
		try(Output out=new Output(prefix+".profile_summary.tsv")){
			out.writer.println("inputs_observed\tcatalog_inputs\tqueried_models\traw_hit_events\tsource_claims\tunavailable_claims\tsource_feature_rows\tretained_rows\texcluded_rows\tevidence_status\tstatus");
			b.clear().append(observed).tab().append(catalog.size()).tab().append(models.size()).tab().append(hits).tab().append(claims).tab();
			b.append(unavailable).tab().append(features).tab().append(retained).tab().append(excluded).tab().append(evidenceState).tab();
			out.writer.print(b.append(observed==catalog.size() ? "PROFILE_COMPLETE_NOT_ADJUDICATED" : "SAMPLE_PROFILE_NOT_ADJUDICATED").nl());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class ModelProfile{
		ModelProfile(final Model model_){
			model=model_;
			bins=new long[(model.length+1)*4];
		}
		long at(final int span){
			final int start=span*4;
			return bins[start]+bins[start+1]+bins[start+2]+bins[start+3];
		}
		/** Nearest-rank weighted quantile; integer decomposition avoids overflow. */
		int quantile(final int percent){
			assert(sum>0 && percent>0 && percent<=100) : "Quantiles require observed events and a positive percentile";
			final long rank=(sum/100)*percent+((sum%100)*percent+99)/100;
			long count=0;
			for(int span=1; span<=model.length; span++){
				count+=at(span);
				if(count>=rank){return span;}
			}
			throw new IllegalStateException("Conserved histogram did not reach its quantile rank");
		}
		/** Smallest span wins a tied modal count; this is a descriptive convention, not a cutoff. */
		int mode(){
			assert(sum>0) : "Zero-hit accessions have no modal span";
			int result=0;
			long maximum=0;
			for(int span=1; span<=model.length; span++){
				final long count=at(span);
				if(count>maximum){
					maximum=count;
					result=span;
				}
			}
			return result;
		}
		final Model model;
		final long[] bins, truncations=new long[4];
		long denominator=-1, sum;
		boolean zeroRow;
	}

	private static final class Availability{
		Availability(final Genome genome_, final ByteBuilder buffer){
			genome=genome_;
			paths=new byte[genome.memberships()][];
			sourceRows=new byte[paths.length][];
			members=new int[paths.length];
			seen=new boolean[paths.length];
			for(int i=0; i<genome.memberships(); i++){
				final byte[] path=BuildNcrnaReference.field(buffer.clear(), genome.membership(i).gff).toBytes();
				int index=0;
				while(index<claims && !Tools.equals(path, paths[index])){index++;}
				if(index==claims){paths[claims++]=path;}
				members[index]++;
			}
		}
		/** Compares serialized UTF-8 bytes; String termEquals is ASCII-only. */
		int find(final LineParser1 parser){
			parser.setBounds(1);
			for(int i=0; i<claims; i++){
				final byte[] expected=paths[i];
				if(expected.length!=parser.b()-parser.a()){continue;}
				int j=0;
				while(j<expected.length && expected[j]==parser.line()[parser.a()+j]){j++;}
				if(j==expected.length){return i;}
			}
			return -1;
		}
		final Genome genome;
		final byte[][] paths, sourceRows;
		final int[] members;
		final boolean[] seen;
		int claims, seenClaims, nonempty, empty, unavailable;
		long features, retained, excluded;
	}

	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private final NcrnaEvidenceCatalog catalog;
	private final Roster roster;
	private final ArrayList<ModelProfile> models=new ArrayList<ModelProfile>();
	private final ObjectMap<String,ModelProfile> byModel=new ObjectMap<String,ModelProfile>(String.class, ModelProfile.class);
	private final ArrayList<Availability> inputs=new ArrayList<Availability>();
	private final ObjectMap<String,Availability> byInput=new ObjectMap<String,Availability>(String.class, Availability.class);
	private int observed;
	private String evidenceState;
	private boolean batchEvidence;
	private long bases, contigs, hits, claims, unavailable, features, retained, excluded, historicalPairs, modelPairs;
	static final String[] INPUT_SUFFIXES={".summary.tsv", ".coverage.tsv", ".sources.tsv.gz"};
	static final String[] OUTPUT_SUFFIXES={".model_profile.tsv", ".input_availability.tsv", ".availability.tsv", ".profile_summary.tsv"};
	private static final String[] CATEGORIES={"NO_AVAILABLE_GFF", "AVAILABLE_ONLY", "MIXED_AVAILABLE_UNAVAILABLE"};
	private static final String[] UNAVAILABLE={"MISSING_GFF", "UNRESOLVED_GFF_PATH", "NOT_REGULAR_GFF", "UNREADABLE_GFF", "UNRESOLVED_GFF_IO"};
}
