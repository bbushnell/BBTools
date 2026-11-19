package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import fileIO.ReadWrite;
import map.IntHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Selected;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Compact consensus/member profiles on saved windows, unpartitioned by default.
 * Follows template.A_SampleByteFile's parser and owner-managed native streams.
 * Reads accepted identity/cohort values; it cannot measure new identities.
 * One member-hashing pass per width supplies every declared consensus weight.
 * Queries use immutable postings and reused primitive score/rank scratch.
 * No genomes, caller, alignment, fitted parameter or resource replacement.
 * @author Brian, Raiden
 */
public final class RankCompactNcrnaWindows{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, RankCompactNcrnaWindows.class, false);
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			final RankDualNcrnaWindows.Config c=new RankDualNcrnaWindows.Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			boolean partitions=false;
			boolean fineWeights=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("panels")){c.panels=b;}
				else if(a.equals("selection")){c.selection=b;}
				else if(a.equals("baseline")){c.baseline=b;}
				else if(a.equals("identities")){c.identities=b;}
				else if(a.equals("members")){c.members=b;}
				else if(a.equals("ledger")){c.ledger=b;}
				else if(a.equals("assignments")){c.assignments=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("summary")){c.summary=b;}
				else if(a.equals("ref")){c.references=b;}
				else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
				else if(a.equals("partitions")){partitions=Parse.parseBoolean(b);}
				else if(a.equals("fineweights")){fineWeights=Parse.parseBoolean(b);}
				else{require(parser.parse(arg, a, b), "Unknown compact saved-window parameter: "+arg);}
			}
			require(parser.out1!=null && parser.out2==null && parser.in1==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Declare complete saved bindings and out=fresh_prefix without genomes, sampling, overwrite or append");
			c.output=parser.out1;
			DynamicDemiLog.setExponent(4);
			require(!partitions || !fineWeights, "Partition geometry and unpartitioned fine-weight studies are separate saved comparisons");
			run(c, partitions, fineWeights);
			pp.outstream.println(fineWeights ? "NCRNA_COMPACT_FINE_WEIGHTS_SAVED_RANKS_COMPLETE_NOT_CALLER_ADOPTED"
				: partitions ? "NCRNA_COMPACT_PARTITIONS_SAVED_RANKS_COMPLETE_NOT_CALLER_ADOPTED"
				: "NCRNA_COMPACT_S1_SAVED_RANKS_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(exponent);
			Shared.closeStream(pp.outstream);
		}
	}

	public static void run(final RankDualNcrnaWindows.Config c){
		run(c, false, false);
	}

	/** The partition pass reuses accepted S1 controls externally; it never
	 * rebuilds those profiles or repeats their scientific saved-window comparison. */
	public static void runPartitions(final RankDualNcrnaWindows.Config c){
		run(c, true, false);
	}

	/** Only previously unmeasured weights at 1,024 bins. Existing 10/100 results
	 * supply those requested rows externally; do not rerun finished comparisons. */
	public static void runFineWeights(final RankDualNcrnaWindows.Config c){
		run(c, false, true);
	}

	private static void run(final RankDualNcrnaWindows.Config c, final boolean partitions, final boolean fineWeights){
		require(c!=null && c.panels!=null && c.selection!=null && c.baseline!=null && c.identities!=null
			&& c.members!=null && c.ledger!=null && c.assignments!=null && c.census!=null
			&& c.summary!=null && c.references!=null && c.output!=null && c.models>0 && c.groups>0,
			"Compact ranks require original cache/member/reference bindings AND accepted identity/cohort files");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Compact saved ranks need assertions, native exponent4 and no blacklist");
		RrnaResourceIO.fresh(c.output, SUFFIXES);
		final int zip=ReadWrite.MAX_ZIP_THREADS();
		final boolean zipSet=ReadWrite.SET_ZIP_THREADS, trim=Shared.TRIM_READ_DESCRIPTION;
		try{
			ReadWrite.setZipThreads(1);
			Shared.TRIM_READ_DESCRIPTION=false;
			new RankCompactNcrnaWindows(c, partitions, fineWeights).process();
		}finally{
			ReadWrite.setZipThreads(zip);
			ReadWrite.SET_ZIP_THREADS=zipSet;
			Shared.TRIM_READ_DESCRIPTION=trim;
		}
	}

	private RankCompactNcrnaWindows(final RankDualNcrnaWindows.Config c_, final boolean partitions_, final boolean fineWeights_){
		c=c_;
		partitions=partitions_;
		fineWeights=fineWeights_;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved Inputs         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			sequences.add(read.bases);
		});
		require(names.size()==c.models, "Original reference order must contain every cached model");
		labels=new byte[c.models][];
		for(int model=0; model<c.models; model++){labels[model]=names.get(model).getBytes(StandardCharsets.UTF_8);}
		final String[] headers=names.toArray(new String[names.size()]);
		consensuses=sequences.toArray(new byte[sequences.size()][]);
		members=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments, c.census, c.summary, headers, consensuses);
		if(!fineWeights){data.readNativeCounts(headers);}
		baseline=RankDualNcrnaWindows.readBaseline(c.baseline, data, c.models);
		cohorts=NcrnaDdlCohorts.read(c.identities, data, labels, consensuses);
		for(final Window window:data.windows){
			require(RankDualNcrnaWindows.invalidBase(window.bases)<0, "Captured windows must retain the accepted ACGTN encoding domain");
		}
		try(Output methods=new Output(c.output+SUFFIXES[0]); Output models=new Output(c.output+SUFFIXES[1]);
				Output loci=new Output(c.output+SUFFIXES[2]); Output summary=new Output(c.output+SUFFIXES[3]);
				Output losses=new Output(c.output+SUFFIXES[4]); Output inputs=new Output(c.output+SUFFIXES[5])){
			final Reports out=new Reports(methods, models, loci, summary, losses, inputs);
			out.inputs.writer.println("panel\tinput_slot\tinput_id\tfasta\ttarget_loci\tnative_common\tcaptured_windows");
			for(final Selected input:data.selected){
				out.inputs.writer.print(b.clear().append(input.panel.nameBytes).tab().append(input.slot).tab().append(input.idBytes)
					.tab().append(input.pathBytes).tab().append(input.parents).tab().append(input.common).tab().append(input.windows).nl());
			}
			reference("native", 0, baseline[0], out);
			reference("historical_hybrid_lowest_key_2048", 2048, baseline[2], out);
			if(fineWeights){
				rankWidth(1024, out, new int[]{2,3,5,8,20,50}, false);
				return;
			}
			for(final int width:new int[]{512, 1024, 2048}){
				if(partitions){
					for(final int shortK:new int[]{0,7,8}){rankPartition(shortK, width, out);}
				}else{rankWidth(width, out);}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Compact Construction    ----------------*/
	/*--------------------------------------------------------------*/

	private void rankWidth(final int width, final Reports out){
		final int[] weights=width==2048 ? new int[]{1000000,0} : new int[]{1000000,0,1,10,100,1000};
		rankWidth(width, out, weights, true);
	}

	/** The fine pass reuses ranking/reporting while omitting completed controls. */
	private void rankWidth(final int width, final Reports out, final int[] weights, final boolean controls){
		final int[] constructionWeights=Arrays.copyOf(weights, weights.length+(controls ? 1 : 0));
		if(controls){constructionWeights[weights.length]=NcrnaDdlProfile.CONSENSUS_WINS;}
		final NcrnaDdlProfile[][] profiles=NcrnaDdlModelProfiles.buildConsensus(members, consensuses, width, constructionWeights);
		final char[][] consensusKeys=controls ? new char[c.models][] : null;
		final NcrnaDdlQueryKeys query=new NcrnaDdlQueryKeys(width, 10);
		final String[] methods=new String[weights.length+(controls ? 2 : 0)];
		for(int variant=0; variant<weights.length; variant++){
			if(weights[variant]==0){methods[variant]="member_highest_"+width;}
			else{methods[variant]="consensus_weight"+weights[variant]+'_'+width;}
		}
		if(controls){
			methods[weights.length]="consensus_only_"+width;
			methods[weights.length+1]="native_first_member_highest_"+width;
		}
		final long[] matches=new long[methods.length];
		final int[] minimum=new int[methods.length], maximum=new int[methods.length];
		Arrays.fill(minimum, Integer.MAX_VALUE);
		long lengthSum=0, collisionSum=0;
		int lengthMin=Integer.MAX_VALUE, lengthMax=0, collisionMin=Integer.MAX_VALUE, collisionMax=0;
		final IntHashSet rawWords=new IntHashSet(256);
		final int[] rawBins=new int[width];
		for(int model=0; model<c.models; model++){
			query.build(consensuses[model], null, false, true);
			final char[] keys=new char[width];
			for(int i=0; i<query.keyCount(); i++){
				final int packed=query.packedKey(i), bucket=packed>>>16;
				keys[bucket]=(char)Math.max(keys[bucket], packed&65535);
			}
			if(controls){consensusKeys[model]=keys;}
			final int occupied=occupied(keys), length=consensuses[model].length;
			final int words=countConsensusWords(consensuses[model], rawWords, rawBins);
			final int collisionBins=countCollisionBins(rawBins);
			int differentBins=0;
			final char[] million=controls ? profiles[0][model].keys() : null;
			final char[] legacy=controls ? profiles[weights.length][model].keys() : null;
			if(controls){
				for(int bucket=0; bucket<width; bucket++){
					if(million[bucket]!=legacy[bucket]){differentBins++;}
				}
			}
			collisionSum+=collisionBins;
			collisionMin=Math.min(collisionMin, collisionBins);
			collisionMax=Math.max(collisionMax, collisionBins);
			lengthSum+=length;
			lengthMin=Math.min(lengthMin, length);
			lengthMax=Math.max(lengthMax, length);
			for(int variant=0; variant<methods.length; variant++){
				final NcrnaDdlProfile profile=variant<weights.length ? profiles[variant][model]
					: variant==weights.length ? null : profiles[1][model];
				final char[] selected=profile==null ? keys : profile.keys();
				final int matching=query.equalBuckets(selected);
				assert(!controls || variant!=0 || members.members(model).length>=1000000 || matching==occupied)
					: "Weight one million exceeds this model's member count, so every occupied consensus bin must hold a consensus key";
				matches[variant]+=matching;
				minimum[variant]=Math.min(minimum[variant], matching);
				maximum[variant]=Math.max(maximum[variant], matching);
				b.clear().append(methods[variant]).tab().append(model).tab().append(labels[model])
					.tab().append(length).tab().append(members.members(model).length).tab().append(occupied).tab().append(matching)
					.tab().append(occupied(selected)).tab().append(2L*width).tab().append(words).tab().append(collisionBins)
					.tab();
				if(controls){b.append(differentBins);}else{b.append("NA");}
				b.tab();
				if(differentBins==0){b.append("NA");}
				else{
					boolean first=true;
					for(int bucket=0; bucket<width; bucket++){
						if(million[bucket]==legacy[bucket]){continue;}
						if(!first){b.append(',');}
						b.append(bucket).append(':').append((int)legacy[bucket]).append('>').append((int)million[bucket]);
						first=false;
					}
				}
				out.models.writer.print(b.nl());
			}
		}
		final DDLIndexSmall[] indices=new DDLIndexSmall[weights.length+(controls ? 1 : 0)];
		for(int variant=0; variant<weights.length; variant++){
			final char[][] keys=new char[c.models][];
			for(int model=0; model<c.models; model++){keys[model]=profiles[variant][model].keys();}
			indices[variant]=new DDLIndexSmall(keys);
		}
		if(controls){indices[weights.length]=new DDLIndexSmall(consensusKeys);}
		final int[][] ranks=new int[methods.length][data.parents.size()];
		for(final int[] row:ranks){Arrays.fill(row, -1);}
		final float[] scores=new float[c.models];
		final int[] counts=new int[c.models];
		final boolean[] acceptable=new boolean[c.models];
		final int[] tuple=new int[4];
		final NcrnaDdlInterleave order=controls ? new NcrnaDdlInterleave(c.models) : null;
		for(final Window window:data.windows){
			query.build(window.bases, null, false, true);
			for(int variant=0; variant<indices.length; variant++){
				query.count(indices[variant], counts);
				for(int model=0; model<c.models; model++){scores[model]=(float)counts[model];}
				if(controls && variant==1){order.build(window.nativeCounts, scores);}
				for(final Edge edge:window.edges){
					if(edge.unique==0){continue;}
					for(int model=0; model<c.models; model++){acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;}
					NcrnaDdlRanks.best(scores, acceptable, tuple);
					final int parent=edge.parent.index;
					ranks[variant][parent]=NcrnaDdlSavedData.minimum(ranks[variant][parent], tuple[0]);
					if(controls && variant==1){
						order.hybridBest(acceptable, tuple);
						ranks[weights.length+1][parent]=NcrnaDdlSavedData.minimum(ranks[weights.length+1][parent], tuple[0]);
					}
				}
			}
		}
		for(int variant=0; variant<methods.length; variant++){
			out.methods.writer.print(b.clear().append(methods[variant]).tab().append(width).tab().append(2L*width)
				.tab().append(2L*width*c.models).tab().append(c.models).tab().append(matches[variant]/(double)c.models, 6)
				.tab().append(minimum[variant]).tab().append(maximum[variant]).tab().append(lengthSum/(double)c.models, 6)
				.tab().append(lengthMin).tab().append(lengthMax).tab().append(collisionSum/(double)c.models, 6)
				.tab().append(collisionMin).tab().append(collisionMax)
				.append("\tUINT16_KEYS_ONLY;POSTINGS_HEADERS_AND_HEAP_EXCLUDED\n"));
			report(methods[variant], width, ranks[variant], out);
		}
	}

	/** Construction diagnostic only: count distinct masked words before key
	 * compression, so two different words with the same encoded key still count.
	 * Full-span enumeration independently rejects N even at the masked position. */
	public static int countConsensusWords(final byte[] sequence, final IntHashSet words, final int[] bins){
		assert(sequence!=null && Integer.bitCount(bins.length)==1) : "Consensus word counts use the selected model and actual power-of-two profile width";
		words.clear();
		Arrays.fill(bins, 0);
		for(int start=0; start+10<=sequence.length; start++){
			int word=0;
			boolean valid=true;
			for(int offset=0; offset<10; offset++){
				final int base=sequence[start+offset]&~32;
				if(base!='A' && base!='C' && base!='G' && base!='T'){
					valid=false;
					break;
				}
				word=(word<<2)|((base>>1)&3);
			}
			word&=~(3<<8);
			if(valid && words.add(word)){
				final long hash=Tools.hash64shift(word^12345L);
				bins[(int)(hash&(bins.length-1))]++;
			}
		}
		return words.size();
	}

	/** Number of bins containing at least two distinct raw masked words. The
	 * same positive-count diagnostic is used by report rows and the fixed fixture. */
	public static int countCollisionBins(final int[] bins){
		assert(bins!=null && bins.length>0) : "Consensus collision statistics require the actual raw-word bin counts";
		int count=0;
		for(final int words:bins){if(words>1){count++;}}
		return count;
	}

	/** One immutable posting index per partition geometry, with no S1 rehash.
	 * Consensus occupancy counts the reserved quarter, not matching member bins. */
	private void rankPartition(final int shortK, final int width, final Reports out){
		final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.buildPartition(members, consensuses, shortK, width);
		final String method=(shortK==0 ? "partition_family_masked10_" : "partition_plain"+shortK+"_masked10_")+width;
		final char[][] keys=new char[c.models][];
		long lengthSum=0, matchingSum=0;
		int lengthMin=Integer.MAX_VALUE, lengthMax=0, matchingMin=Integer.MAX_VALUE, matchingMax=0;
		for(int model=0; model<c.models; model++){
			keys[model]=profiles[model].keys();
			final int length=consensuses[model].length;
			int consensusBins=0;
			for(int bucket=3*(width/4); bucket<width; bucket++){
				if(keys[model][bucket]!=0){consensusBins++;}
				//finishPartition rejects member keys here before installing consensus maxima.
			}
			lengthSum+=length;
			lengthMin=Math.min(lengthMin, length);
			lengthMax=Math.max(lengthMax, length);
			matchingSum+=consensusBins;
			matchingMin=Math.min(matchingMin, consensusBins);
			matchingMax=Math.max(matchingMax, consensusBins);
			out.models.writer.print(b.clear().append(method).tab().append(model).tab().append(labels[model])
				.tab().append(length).tab().append(profiles[model].members).tab().append(consensusBins).tab().append(consensusBins)
				.tab().append(occupied(keys[model])).tab().append(2L*width).append("\tNA\tNA\tNA\tNA\n"));
		}
		final DDLIndexSmall index=new DDLIndexSmall(keys);
		final NcrnaDdlPartitionKeys query=new NcrnaDdlPartitionKeys(shortK, width);
		final int[] ranks=new int[data.parents.size()], tuple=new int[4];
		Arrays.fill(ranks, -1);
		final float[] scores=new float[c.models];
		final int[] counts=new int[c.models];
		final boolean[] acceptable=new boolean[c.models];
		for(final Window window:data.windows){
			query.build(window.bases);
			query.count(index, counts);
			for(int model=0; model<c.models; model++){scores[model]=(float)counts[model];}
			for(final Edge edge:window.edges){
				if(edge.unique==0){continue;}
				for(int model=0; model<c.models; model++){acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;}
				NcrnaDdlRanks.best(scores, acceptable, tuple);
				final int parent=edge.parent.index;
				ranks[parent]=NcrnaDdlSavedData.minimum(ranks[parent], tuple[0]);
			}
		}
		out.methods.writer.print(b.clear().append(method).tab().append(width).tab().append(2L*width)
			.tab().append(2L*width*c.models).tab().append(c.models).tab().append(matchingSum/(double)c.models, 6)
			.tab().append(matchingMin).tab().append(matchingMax).tab().append(lengthSum/(double)c.models, 6)
			.tab().append(lengthMin).tab().append(lengthMax)
			.append("\tNA\tNA\tNA\tUINT16_KEYS_ONLY;CONSENSUS_STATS_RESERVED_QUARTER;POSTINGS_HEADERS_AND_HEAP_EXCLUDED\n"));
		report(method, width, ranks, out);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reporting          ----------------*/
	/*--------------------------------------------------------------*/

	private void reference(final String method, final int width, final int[] ranks, final Reports out){
		assert(width==0 || width==2048) : "Only the frozen native and historical hybrid controls have unavailable construction statistics";
		out.methods.writer.print(b.clear().append(method).tab().append(width).tab().append(2L*width)
			.tab().append(2L*width*c.models).tab().append(c.models)
			.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNATIVE_KEYS_ARE_IMPLICIT_DENSE_ADDRESSES;DDL_KEYS_ONLY;POSTINGS_AND_HEAP_EXCLUDED\n"));
		report(method, width, ranks, out);
	}

	private void report(final String method, final int width, final int[] ranks, final Reports out){
		assert(ranks.length==data.parents.size()) : "Every method must retain the original full target roster";
		for(final Parent parent:data.parents){
			out.loci.writer.print(b.clear().append(method).tab().append(parent.input.panel.nameBytes).tab().append(parent.input.slot)
				.tab().append(parent.input.idBytes).tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum])
				.tab().append(NcrnaDdlCohorts.BINS[cohorts.bins[parent.index]]).tab().append(parent.nativeRank>0)
				.tab().append(ranks[parent.index]).tab().append(parent.reason).nl());
		}
		final int[][] reachable=new int[2][data.selected.size()];
		final int[] after=new int[data.selected.size()];
		for(int state=0; state<3; state++){
			for(int bin=0; bin<NcrnaDdlCohorts.BINS.length; bin++){
				for(final int[] row:reachable){Arrays.fill(row, 0);}
				Arrays.fill(after, 0);
				int targets=0, common=0, within10=0, within50=0;
				long sum=0;
				for(final Parent parent:data.parents){
					final int p=parent.index;
					if(parent.stratum!=state || bin<5 && cohorts.bins[p]!=bin){continue;}
					targets++;
					if(parent.nativeRank<=0){continue;}
					common++;
					final int rank=ranks[p];
					require(rank>0 && rank<=c.models, "Every compact method must rank every native-common parent");
					sum+=rank;
					if(rank<=10){within10++;}
					if(rank<=50){
						within50++;
						after[parent.input.index]++;
					}
					for(int reference=0; reference<2; reference++){
						final int prior=baseline[reference==0 ? 0 : 2][p];
						require(prior>0, "The comparison must retain both accepted controls for every native-common parent");
						if(prior<=50){reachable[reference][parent.input.index]++;}
					}
				}
				b.clear().append(method).tab().append(NcrnaDdlSavedData.STRATA[state]).tab().append(NcrnaDdlCohorts.BINS[bin])
					.tab().append(targets).tab().append(common).tab().append(sum).tab().append(within10).tab().append(within50)
					.tab().append(2L*width).tab().append(2L*width*c.models);
				for(int reference=0; reference<2; reference++){
					int eligible=0, lost=0;
					for(final Selected input:data.selected){
						if(reachable[reference][input.index]>0){
							eligible++;
							if(after[input.index]==0){lost++;}
						}
					}
					b.tab().append(eligible).tab().append(lost);
				}
				out.summary.writer.print(b.append("\tSAVED_NATIVE_COMMON_LOCI_NOT_CALLER_DETECTION\n"));
				for(int reference=0; reference<2; reference++){
					for(final Selected input:data.selected){
						if(reachable[reference][input.index]==0 || after[input.index]>0){continue;}
						out.losses.writer.print(b.clear().append(method).tab().append(NcrnaDdlSavedData.STRATA[state])
							.tab().append(NcrnaDdlCohorts.BINS[bin]).tab().append(reference==0 ? "native" : "historical_hybrid_lowest_key_2048")
							.tab().append(input.panel.nameBytes).tab().append(input.slot).tab().append(input.idBytes).tab().append(input.pathBytes)
							.tab().append(reachable[reference][input.index]).nl());
					}
				}
			}
		}
	}

	private static int occupied(final char[] keys){
		assert(keys!=null) : "Profile occupancy requires the actual selected key payload";
		int count=0;
		for(final char key:keys){if(key!=0){count++;}}
		return count;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/** Borrowed writers; the owner's try-with-resources handles every close. */
	private static final class Reports{
		Reports(final Output methods_, final Output models_, final Output loci_, final Output summary_,
				final Output losses_, final Output inputs_){
			methods=methods_;
			models=models_;
			loci=loci_;
			summary=summary_;
			losses=losses_;
			inputs=inputs_;
			methods.writer.println("method\ttotal_bins\tkey_bytes_per_model\tkey_bytes_all_models\tmodels\tconsensus_matching_bins_mean\tconsensus_matching_bins_min\tconsensus_matching_bins_max\tconsensus_length_mean\tconsensus_length_min\tconsensus_length_max\tconsensus_collision_bins_mean\tconsensus_collision_bins_min\tconsensus_collision_bins_max\tstorage_scope");
			models.writer.println("method\tmodel_index\tmodel\tconsensus_length\tactual_members\tconsensus_occupied_bins\tconsensus_matching_selected_bins\tselected_occupied_bins\tkey_bytes\tconsensus_distinct_masked10_words\tconsensus_collision_bins\tmillion_vs_always_different_bins\tdifferent_bins_bucket_old_new");
			loci.writer.println("method\tpanel\tinput_slot\tinput_id\tparent\tseen_stratum\tidentity_bin\tnative_common\tbest_rank\tsource_disposition");
			summary.writer.println("method\tseen_stratum\tidentity_bin\ttarget_loci\tcommon_loci\trank_sum\twithin10\twithin50\tkey_bytes_per_model\tkey_bytes_all_models\tnative_reachable_inputs\tnative_last_losses_at50\thistorical_hybrid_reachable_inputs\thistorical_hybrid_last_losses_at50\tinterpretation");
			losses.writer.println("method\tseen_stratum\tidentity_bin\treference\tpanel\tinput_slot\tinput_id\tfasta\treference_reachable_loci");
		}
		final Output methods, models, loci, summary, losses, inputs;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final RankDualNcrnaWindows.Config c;
	private final boolean partitions, fineWeights;
	private final ByteBuilder b=new ByteBuilder();
	private NcrnaDdlSavedData data;
	private NcrnaDdlMembers members;
	private NcrnaDdlCohorts cohorts;
	private byte[][] labels, consensuses;
	private int[][] baseline;
	private static final String[] SUFFIXES={".methods.tsv", ".profile_models.tsv", ".loci.tsv.gz", ".summary.tsv", ".last_losses.tsv", ".inputs.tsv"};
}
