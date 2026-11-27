package prok;

import java.util.ArrayList;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Builds source-bound cobalamin member-voted profiles once, before callers.
 * Forked from template/A_SampleByteFile's argument/validation/process lifecycle;
 * existing native member/FASTA owners handle the joined construction inputs.
 * Voting is one owner pass with reused sketch scratch, not per-caller work.
 * No genome, alignment, verification or caller invocation.
 * @author Brian, Raiden, Ganyu
 */
public final class BuildNcrnaDdlProfiles{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BuildNcrnaDdlProfiles.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			if(appendRequested(pp.args)){
				appendProfiles(pp.args);
				pp.outstream.println("NCRNA_DDL_PROFILE_APPEND_PASS_NO_VOTING");
				return;
			}
			final Config c=parse(pp.args);
			DynamicDemiLog.setExponent(4);
			build(c);
			timer.stop();
			pp.outstream.println("Profile models: "+c.models+"; kept member votes: "+c.members);
			if(c.refit){pp.outstream.println("Profile fit: all eligible originals; no panel, alias or duplicate-voter exclusions");}
			pp.outstream.println("Time: "+timer);
			pp.outstream.println(c.prepareOnly ? "NCRNA_DDL_TRAINING_SAMPLE_COMPLETE" : "NCRNA_DDL_PROFILE_RESOURCE_BUILD_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Explicit resource append rejects construction options rather than silently
	 * applying a different voting population. The ordinary CLI path is unchanged. */
	private static boolean appendRequested(final String[] args){
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair[0].equalsIgnoreCase("appendprofiles")){return true;}
		}
		return false;
	}

	private static void appendProfiles(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		String firstRef=null, secondRef=null;
		int[] expected=null;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("appendprofiles")){require(Parse.parseBoolean(b), "Use appendprofiles=t for explicit profile append");}
			else if(a.equals("ref")){firstRef=b;}
			else if(a.equals("ref2")){secondRef=b;}
			else if(a.equals("expectedmodels")){expected=Parse.parseIntArray(b, ",");}
			else{require(parser.parse(arg, a, b), "Unknown profile-append option (no construction flags): "+arg);}
		}
		require(parser.in1!=null && parser.in2!=null && parser.out1!=null && parser.out2==null
			&& firstRef!=null && secondRef!=null && expected!=null && expected.length==2 && expected[0]>0 && expected[1]>0
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Profile append requires in=old.tsv in2=new.tsv ref=old.fa ref2=new.fa expectedmodels=OLD,NEW out=fresh.tsv without sampling");
		require(Tools.testInputFiles(false, true, parser.in1, parser.in2, firstRef, secondRef), "All profile and consensus inputs must be readable and distinct");
		RrnaResourceIO.fresh(parser.out1, "");
		Shared.setThreads(1);
		final boolean trim=Shared.TRIM_READ_DESCRIPTION;
		try{
			Shared.TRIM_READ_DESCRIPTION=false;
			final NcrnaSupport.Library first=NcrnaSupport.loadLibrary(firstRef);
			final NcrnaSupport.Library second=NcrnaSupport.loadLibrary(secondRef);
			require(first.names.length==expected[0] && second.names.length==expected[1], "Both consensus model counts must match their declared append axes");
			NcrnaDdlProfileResource.appendOrdering(parser.in1, parser.in2, parser.out1,
				first.names, first.bases, second.names, second.bases);
		}finally{Shared.TRIM_READ_DESCRIPTION=trim;}
	}

	private static Config parse(final String[] args){
		assert(args!=null) : "PreParser must supply the expanded native argument list";
		final Config c=new Config();
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("ref")){c.references=b;}
			else if(a.equals("ledger")){c.ledger=b;}
			else if(a.equals("assignments")){c.assignments=b;}
			else if(a.equals("additions")){c.additions=b;}
			else if(a.equals("census")){c.census=b;}
			else if(a.equals("summary")){c.summary=b;}
			else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmembers")){
				c.members=Parse.parseIntKMG(b);
				c.membersSet=true;
			}
			else if(a.equals("dual")){c.dual=Parse.parseBoolean(b);}
			else if(a.equals("partitionbins")){c.partitionBins=Parse.parseIntKMG(b);}
			else if(a.equals("selected")){c.selected=Parse.parseBoolean(b);}
			else if(a.equals("orderonly")){c.orderOnly=Parse.parseBoolean(b);}
			else if(a.equals("refit")){c.refit=Parse.parseBoolean(b);}
			else if(a.equals("suppressconsensus")){c.suppressConsensus=Parse.parseBoolean(b);}
			else if(a.equals("consensusbonus")){c.consensusBonus=Parse.parseIntKMG(b);}
			else if(a.equals("consensusrule")){c.consensusRule=NcrnaDdlProfileResource.parseBonusRule(b);}
			else if(a.equals("prepareonly")){c.prepareOnly=Parse.parseBoolean(b);}
			else if(a.equals("seed")){c.seeds=new long[]{Parse.parseKMG(b)};}
			else if(a.equals("seeds")){c.seeds=Parse.parseLongArray(b);}
			else if(a.equals("manifest")){c.manifest=b;}
			else if(a.equals("panel")){c.panel=b;}
			else if(a.equals("sources")){c.sources=b;}
			else if(a.equals("protectedpriors")){c.protectedPriors=b;}
			else if(a.equals("expectedpriors")){c.expectedPriors=Parse.parseIntKMG(b);}
			else if(a.equals("aliases")){c.aliases=b;}
			else if(a.equals("trainsample")){c.trainSample=Parse.parseIntKMG(b);}
			else if(a.equals("expectedpanel")){c.panelInputs=Parse.parseIntKMG(b);}
			else if(a.equals("sampleseed")){c.sampleSeed=Parse.parseKMG(b);}
			else if(a.equals("expectedinput")){c.inputMembers=Parse.parseIntKMG(b);}
			else if(a.equals("paddingblacklist")){c.paddingBlacklist=b;}
			else if(a.equals("paddingfit")){c.paddingFit=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown profile-build parameter: "+arg);}
		}
		require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Use in=construction_members.fa out=fresh_profiles.tsv with no sampling, overwrite or append");
		c.fasta=parser.in1;
		c.output=parser.out1;
		return c;
	}

	private static void validate(final Config c){
		require(c!=null && c.fasta!=null && c.references!=null && c.assignments!=null
			&& c.output!=null && c.models>0 && (c.members>0 || c.orderOnly && c.members==0),
			"Declare the complete frozen construction join, reference and positive expected model/member counts");
		final boolean roster=c.selected || c.orderOnly;
		require(!c.refit || c.orderOnly && !c.prepareOnly && c.panel==null,
			"refit=t requires orderonly=t, omits panel= and cannot prepare a held-out training sample");
		require(c.protectedPriors==null ? c.expectedPriors==0 : c.orderOnly && !c.prepareOnly
			&& c.expectedPriors>0 && c.inputMembers>c.expectedPriors && Tools.testInputFiles(false, true, c.protectedPriors),
			"Protected priors require explicit ordering-only construction, companion/count and the complete augmented input population");
		NcrnaDdlProfileResource.validateOrderingBonus(c.consensusBonus, c.consensusRule);
		require(c.consensusBonus==0 || c.orderOnly && !c.prepareOnly && !c.suppressConsensus,
			"An in-profile consensus bonus is an explicit ordering-only build, never suppression, preparation or a legacy format");
		require(!c.suppressConsensus || c.orderOnly && !c.prepareOnly && c.partitionBins==1024,
			"suppressconsensus=t requires actual ordering-only 1024-bin profile construction");
		require(c.additions==null || c.orderOnly && !c.prepareOnly && Tools.testInputFiles(false, true, c.additions),
			"Explicit singleton additions require ordering-only profiles and their readable producer ledger");
		require(!c.prepareOnly || c.orderOnly, "prepareonly=t is the explicit current-roster holdout/sample step");
		require(c.seeds!=null && c.seeds.length>0, "Declare at least one explicit sketch seed");
		final LongHashSet uniqueSeeds=new LongHashSet(c.seeds.length*2);
		for(final long seed:c.seeds){require(seed>=0 && uniqueSeeds.add(seed), "Sketch seeds must be explicit, nonnegative and distinct");}
		require(c.orderOnly || c.seeds.length==1 && c.seeds[0]==12345,
			"Legacy profile formats retain seed12345; seed variation requires orderonly=t");
		if(roster){
			require(!Shared.AMINO_IN && !Read.IUPAC_TO_N && !Read.U_TO_T && !Read.DOT_DASH_X_TO_N && !Read.LOWER_CASE_TO_N,
				"Selected roster voting preserves actual nucleotide IUPAC content without lossy input conversion");
			require(c.inputMembers>0 && (c.orderOnly || c.inputMembers>=c.members) && c.ledger==null && c.census==null && c.summary==null
				&& (c.prepareOnly || c.paddingBlacklist!=null)
				&& (c.partitionBins==512 || c.partitionBins==1024) && !c.dual,
				"selected=t uses the new roster FASTA/assignments only, explicit expectedinput and partitionbins512/1024, and explicit paddingblacklist=none or a scoped table");
			if(!c.prepareOnly){
				if("none".equals(c.paddingBlacklist)){require(c.paddingFit==null, "An empty padding set has no learned fit scope");}
				else{require(c.paddingFit!=null && !c.paddingFit.isEmpty() && Tools.testInputFiles(false, true, c.paddingBlacklist), "A learned padding set requires its actual table and explicit scope");}
			}
		}else{
			require(c.ledger!=null && c.census!=null && c.summary!=null && c.paddingBlacklist==null && c.paddingFit==null && c.inputMembers<0,
				"Legacy layouts retain their original full construction join; roster-only flags require selected=t");
		}
		if(c.orderOnly){
			require(!c.selected && c.manifest!=null && (c.refit || c.panel!=null) && c.sources!=null && c.aliases!=null
				&& Tools.testInputFiles(false, true, c.manifest, c.panel, c.sources, c.aliases),
				"Ordering-only profiles require explicit formal manifest, evaluation panel, source ledger and reviewed alias map");
			require(c.refit || c.trainSample>0 && c.sampleSeed>=0 && c.panelInputs>0,
				"Training sampling needs explicit positive sample/panel sizes and a nonnegative seed");
			if(c.prepareOnly){
				require(c.paddingBlacklist==null && c.paddingFit==null
					&& Tools.testOutputFiles(false, false, false, c.output+".training.tsv", c.output+".holdout.tsv", c.output+".holdout.tsv.genomes.tsv"),
					"The preparation step writes fresh exclusions/sample only, before padding exists");
			}else{
				for(final long seed:c.seeds){
					require(Tools.testOutputFiles(false, false, false, profilePath(c, seed)), "Every seeded profile output must be fresh");
				}
				require(Tools.testOutputFiles(false, false, false, c.output+".members.tsv"), "The model-voter audit must be fresh");
			}
		}else{
			require(c.manifest==null && c.panel==null && c.sources==null && c.aliases==null,
				"Current holdout metadata requires orderonly=t rather than silently changing a legacy population");
		}
		require((c.partitionBins==0 || c.partitionBins==512 || c.partitionBins==1024)
			&& !(c.dual && c.partitionBins!=0),
			"Use partitionbins=512or1024 for the selected compact design; dual=t is a separate historical geometry");
		require(Tools.testInputFiles(false, true, c.fasta, c.references, c.ledger, c.assignments, c.census, c.summary),
			"Every frozen construction/reference input must be readable");
		require(Tools.testOutputFiles(false, false, false, c.output), "Profile output must be fresh and writable");
		require(Tools.testForDuplicateFiles(true, c.fasta, c.references, c.ledger, c.assignments, c.census, c.summary, c.output),
			"Construction/reference inputs and profile output must have distinct paths");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Profile construction requires assertions, fixed exponent4 and no DDL blacklist");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	public static void build(final Config c){
		validate(c);
		final boolean previousTrim=Shared.TRIM_READ_DESCRIPTION;
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		try{
			Shared.TRIM_READ_DESCRIPTION=false;
			NcrnaDdlSavedData.sequences(c.references, read->{
				names.add(read.id);
				sequences.add(read.bases);
			});
			require(names.size()==c.models, "The complete actual reference roster must match expectedmodels");
			final String[] ids=names.toArray(new String[names.size()]);
			final byte[][] bases=sequences.toArray(new byte[sequences.size()][]);
			final String family=c.orderOnly ? NcrnaDdlMembers.memberFamily(ids) : null;
			final String audit=c.output+(c.refit ? ".refit.tsv" : ".holdout.tsv");
			if(c.orderOnly && !c.prepareOnly && (c.refit || !family.equals("RF00174"))){
				require(Tools.testOutputFiles(false, false, false, audit,
					audit+".genomes.tsv", c.output+".dispositions.tsv"),
					"All minor holdout/disposition outputs must be fresh before any seeded profile is written");
			}
			final NcrnaDdlHoldout holdout=c.orderOnly ? new NcrnaDdlHoldout(c.manifest, c.panel, c.sources,
				c.aliases, family, c.protectedPriors, c.expectedPriors, c.refit) : null;
			if(holdout!=null && !c.refit){require(holdout.panelGenomes()==c.panelInputs, "The excluded panel must match expectedpanel exactly");}
			final NcrnaDdlMembers members=c.orderOnly
				? NcrnaDdlMembers.loadSubsetRoster(c.fasta, c.assignments, c.additions, ids, bases, c.inputMembers, holdout)
				: c.selected
				? NcrnaDdlMembers.loadRoster(c.fasta, c.assignments, ids, bases, c.inputMembers, c.members)
				: NcrnaDdlMembers.load(c.fasta, c.ledger, c.assignments, c.census, c.summary, ids, bases);
			require(c.orderOnly && !c.membersSet || members.votingMembers()==c.members,
				"The retained model-vote incidence population must match any declared expectedmembers");
			c.members=members.votingMembers();
			if(c.prepareOnly){
				holdout.writeSample(c.output+".training.tsv", c.trainSample, c.sampleSeed);
				holdout.writeAudit(c.output+".holdout.tsv");
				return;
			}
			if(c.orderOnly){
				final LongHashSet blacklist="none".equals(c.paddingBlacklist) ? new LongHashSet(4)
					: NcrnaDdlSketchBuilder.readBlacklist(c.paddingBlacklist, c.paddingFit);
				for(final long seed:c.seeds){buildOrdering(c, members, ids, bases, blacklist, seed);}
				writeMembers(c.output+".members.tsv", members, ids);
				if(c.refit || !family.equals("RF00174")){
					holdout.writeAudit(audit);
					writeDispositions(c.output+".dispositions.tsv", members);
				}
			}else if(c.selected){
				final LongHashSet blacklist="none".equals(c.paddingBlacklist) ? new LongHashSet(4)
					: NcrnaDdlSketchBuilder.readBlacklist(c.paddingBlacklist, c.paddingFit);
				final NcrnaDdlSketchBuilder builder=new NcrnaDdlSketchBuilder(c.partitionBins, 12345, blacklist);
				final DynamicDemiLog voter=builder.newSketch();
				final NcrnaDdlSketchBuilder.Prefix prefix=new NcrnaDdlSketchBuilder.Prefix(c.partitionBins);
				final NcrnaDdlLayoutVotes votes=new NcrnaDdlLayoutVotes(c.partitionBins);
				final char[][] keys=new char[c.models][];
				final int[] counts=new int[c.models];
				for(int model=0; model<c.models; model++){
					votes.reset();
					for(final byte[] member:members.members(model)){
						builder.fill(member, voter);
						votes.add(voter.maxArray());
					}
					builder.fillPrefix(bases[model], prefix);
					final NcrnaDdlLayoutVotes.Result result=votes.finish(prefix, 100, true);
					keys[model]=result.keys;
					counts[model]=result.members;
				}
				NcrnaDdlProfileResource.writeSelected(c.output, keys, counts, ids, bases, c.partitionBins, blacklist);
			}else if(c.partitionBins!=0){
				final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.buildPartition(members, bases, 7, c.partitionBins);
				NcrnaDdlProfileResource.writeCompact(c.output, profiles, ids, bases, c.partitionBins);
			}else{
				final NcrnaDdlProfile[] profiles=c.dual
					? NcrnaDdlModelProfiles.buildDual(members, 9, 10, NcrnaDdlProfileResource.DUAL_BUCKETS, true)
					: NcrnaDdlModelProfiles.build(members, NcrnaDdlProfileResource.K, NcrnaDdlProfileResource.BUCKETS, false, true);
				NcrnaDdlProfileResource.write(c.output, profiles, ids, bases, c.dual);
			}
		}finally{Shared.TRIM_READ_DESCRIPTION=previousTrim;}
	}

	private static String profilePath(final Config c, final long seed){
		return c.seeds.length==1 ? c.output : c.output+"."+seed+".tsv";
	}

	private static void buildOrdering(final Config c, final NcrnaDdlMembers members, final String[] ids,
			final byte[][] bases, final LongHashSet blacklist, final long seed){
		assert(c.orderOnly && members.models()==ids.length) : "Ordering profiles keep the fixed current member-model axis";
		final NcrnaDdlSketchBuilder builder=new NcrnaDdlSketchBuilder(c.partitionBins, seed, blacklist);
		final DynamicDemiLog voter=builder.newSketch();
		final NcrnaDdlLayoutVotes votes=new NcrnaDdlLayoutVotes(c.partitionBins);
		final char[][] keys=new char[ids.length][];
		final int[] counts=new int[ids.length];
		final char[] consensus=c.suppressConsensus || c.consensusRule==NcrnaDdlProfileResource.BONUS_MAXIMA
			? new char[c.partitionBins] : null;
		final NcrnaDdlSketchBuilder.Prefix prefix=c.consensusBonus>0 ? new NcrnaDdlSketchBuilder.Prefix(c.partitionBins) : null;
		for(int model=0; model<ids.length; model++){
			votes.reset();
			if(consensus!=null){builder.fillKeys(bases[model], consensus);}
			for(final byte[] member:members.members(model)){
				builder.fill(member, voter);
				if(!c.suppressConsensus){votes.add(voter.maxArray());}
				else{votes.addSuppressed(voter.maxArray(), consensus);}
			}
			if(prefix!=null){
				if(c.consensusRule==NcrnaDdlProfileResource.BONUS_MAXIMA){prefix.setMaxima(consensus);}
				else{builder.fillPrefix(bases[model], prefix);}
			}
			final NcrnaDdlLayoutVotes.Result result=prefix==null ? votes.finish() : votes.finishConsensusBonus(prefix, c.consensusBonus);
			keys[model]=result.keys;
			counts[model]=result.members;
		}
		if(c.consensusBonus>0){
			NcrnaDdlProfileResource.writeBonusOrdering(profilePath(c, seed), keys, counts, ids, bases,
				c.partitionBins, blacklist, seed, c.consensusBonus, c.consensusRule);
		}else if(c.suppressConsensus){
			NcrnaDdlProfileResource.writeSuppressedOrdering(profilePath(c, seed), keys, counts, ids, bases, c.partitionBins, blacklist, seed);
		}else{
			NcrnaDdlProfileResource.writeOrdering(profilePath(c, seed), keys, counts, ids, bases, c.partitionBins, blacklist, seed);
		}
	}

	private static void writeMembers(final String path, final NcrnaDdlMembers members, final String[] names){
		RrnaResourceIO.fresh(path, "");
		try(Output out=new Output(path)){
			out.writer.println("model_index\tmodel\tassigned_before_holdout\tkept_voters\texcluded_voters");
			final ByteBuilder row=new ByteBuilder();
			for(int model=0; model<names.length; model++){
				final int kept=members.members(model).length, assigned=members.assignedMembers(model);
				out.writer.print(row.clear().append(model).tab().append(RrnaResourceIO.first(names[model]))
					.tab().append(assigned).tab().append(kept).tab().append(assigned-kept).nl());
			}
		}
	}

	/** Unique original dispositions and additional vote edges have different denominators. */
	private static void writeDispositions(final String path, final NcrnaDdlMembers members){
		RrnaResourceIO.fresh(path, "");
		try(Output out=new Output(path)){
			out.writer.println("original_genes\tretained_original_genes\texcluded_original_genes\tremoved_cluster_genes\tsmall_genes\torphan_genes\tadded_vote_edges\texcluded_added_edges\tretained_vote_edges");
			out.writer.print(new ByteBuilder().append(members.inputMembers()).tab().append(members.keptMembers())
				.tab().append(members.excludedMembers()).tab().append(members.removedMembers())
				.tab().append(members.smallMembers()).tab().append(members.orphanMembers())
				.tab().append(members.addedMembers()).tab().append(members.excludedAddedMembers())
				.tab().append(members.votingMembers()).nl());
		}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String fasta, references, ledger, assignments, additions, census, summary, output, paddingBlacklist, paddingFit;
		public String manifest, panel, sources, aliases, protectedPriors;
		public int expectedPriors;
		public int models=508, members=85726, partitionBins, inputMembers=-1;
		public int trainSample=4000, panelInputs=1000;
		public int consensusBonus, consensusRule=NcrnaDdlProfileResource.BONUS_NONE;
		public long sampleSeed=20261119;
		public long[] seeds={12345};
		public boolean dual, selected, orderOnly, prepareOnly, membersSet, suppressConsensus, refit;
	}
}
