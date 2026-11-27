package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CallGenes;
import prok.CobalaminCallerResources;
import prok.CobalaminCmAcceptance;
import prok.CobalaminCmEndpointRefiner;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.GeneModel;
import prok.NcrnaClanScavenger;
import prok.NcrnaCpuProfile;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaMatchedBaseWinner;
import prok.NcrnaOverlapWindowBuilder;
import prok.NcrnaPathScoreConfig;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import shared.Shared;

/** Actual serving-group load and ordinary lazy-worker setup versus the explicit
 * development setup on bounded synthetic consensus probes. No genome benchmark.
 * Negative packages contain data symlinks only, never a program/class copy.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCallerResourcesTest{
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminCallerResourcesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String scoreCell="padding";
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
				if(a.equals("scorecell")){scoreCell=b;}
				else{check(parser.parse(arg, a, b), "Use in=serving_group out=fresh_fixture_directory scorecell=padding|greedy");}
			}
			check("padding".equals(scoreCell) || "greedy".equals(scoreCell), "The finite comparison binds its explicit family-only or saved-greedy cell");
			check(parser.in1!=null && parser.out1!=null && !parser.overwrite && !parser.append,
				"The fixture reads one complete serving group and owns a fresh data-only output directory");
			check(!NcrnaCpuProfile.enabled(), "The ordinary-entry fixture must start without CPU profiling");
			final File root=new File(parser.out1), source=new File(parser.in1).getCanonicalFile();
			check(root.mkdir(), "The fixture directory must be fresh");
			final ArrayList<NcrnaFamily> previous=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
			final CobalaminCallerResources old=GeneCaller.cobalaminCallerResources;
			final boolean enabled=CallGenes.COBALAMIN_ENABLED;
			final boolean ncrna=CallGenes.NCRNA_FAMILIES_ENABLED;
			final String oldRoot=CallGenes.COBALAMIN_RESOURCE_ROOT;
			try{
				final NcrnaFamilyLoader.Options disabled=new NcrnaFamilyLoader.Options();
				disabled.enabled=false;
				disabled.resources.root=new File(root, "absent_disabled_group").getPath();
				check(NcrnaFamilyLoader.load(FAMILIES[0], disabled, new FamilyResourceResolver())==null,
					"A disabled cobalamin selection must not resolve the new policy or any payload");
				load(source);
				check(GeneCaller.cobalaminCallerResources!=null && GeneCaller.ncrnaFamilies.size()==3,
					"Ordinary CallGenes loading must publish the complete packaged group and serving configuration");
				final CobalaminCallerResources loaded=GeneCaller.cobalaminCallerResources;
				final NcrnaFamily firstMember=GeneCaller.ncrnaFamilies.get(0);
				GeneCaller.ncrnaFamilies.clear();
				GeneCaller.ncrnaFamilies.add(firstMember);
				caller().ncrnaScavenger(0);
				check(GeneCaller.cobalaminCallerResources==null, "An explicitly reset and replaced roster cannot retain the previous group's serving cache");
				GeneCaller.ncrnaFamilies.clear();
				GeneCaller.cobalaminCallerResources=loaded;
				final File missing=cloneGroup(source, new File(root, "reset_missing_policy"), "missing_policy");
				boolean staleRejected=false;
				try{load(missing);}catch(IllegalArgumentException expected){staleRejected=true;}
				check(staleRejected && GeneCaller.cobalaminCallerResources==null,
					"A failed reload after a roster reset must discard stale serving resources before resolving new assets");
				resetOptions(source, root);
				load(source);
				final NcrnaFamily[] families=GeneCaller.ncrnaFamilies.toArray(new NcrnaFamily[3]);
				final NcrnaPathScoreConfig.Cell scores=bindScores(source, families, scoreCell);
				final byte[] probe=probe(families[0], scores, scoreCell.equals("greedy"));
				for(int strand=0; strand<2; strand++){
					final GeneCaller ordinary=caller();
					final NcrnaClanScavenger actual=(NcrnaClanScavenger)ordinary.ncrnaScavenger(0);
					check(!actual.hasNativeRankingIndex() && actual.modelPatience()==4 && actual.nativeRescue()==null,
						"The ordinary lazy worker must select the merged index-free patient walk");
					final NcrnaClanScavenger control=control(source, families, scores);
					final ArrayList<Orf> first=actual.scavenge("synthetic_consensus", probe.clone(), strand, new ArrayList<int[]>());
					final ArrayList<Orf> second=control.scavenge("synthetic_consensus", probe.clone(), strand, new ArrayList<int[]>());
					check(!first.isEmpty() && first.size()==second.size(), "Both setups must retain the positive synthetic consensus probe");
					boolean changedWinner=!scoreCell.equals("greedy");
					for(int i=0; i<first.size(); i++){
						final Orf a=first.get(i), b=second.get(i);
						check(a.start==b.start && a.stop==b.stop && a.strand==b.strand && a.trnaModel.equals(b.trnaModel)
							&& Float.floatToRawIntBits(a.orfScore)==Float.floatToRawIntBits(b.orfScore),
							"Ordinary setup must preserve explicit-driver member/model/endpoints/path-score bits on both strands");
						for(int model=0; model<families[0].modelNames.length; model++){
							if(families[0].modelNames[model].equals(a.trnaModel) && scores.multiplier(0, model)!=2f){changedWinner=true;}
						}
					}
					check(changedWinner, "Greedy score parity must exercise an emitted model with a changed coefficient, not an unaffected default");
				}
				for(final String fault:new String[]{"missing_profiles", "missing_policy", "wrong_profiles"}){
					final File broken=cloneGroup(source, new File(root, fault), fault);
					boolean rejected=false;
					String reason="No exception was thrown";
					try{load(broken);}catch(IllegalArgumentException expected){
						reason=expected.toString();
						final String needle=fault.equals("missing_policy") ? "caller_policy" : fault.equals("missing_profiles") ? "ddl_profiles"
							: "Ordering-only profiles require the declared member";
						for(Throwable cause=expected; cause!=null; cause=cause.getCause()){
							rejected|=cause.getMessage()!=null && cause.getMessage().contains(needle);
						}
					}
					check(rejected && GeneCaller.ncrnaFamilies.isEmpty() && GeneCaller.cobalaminCallerResources==null,
						"An incomplete or mismatched new group must fail before publishing any member: "+fault+"; "+reason);
				}
				final File withoutOldTips=cloneGroup(source, new File(root, "unused_main_tips"), "unused_main_tips");
				load(withoutOldTips);
				check(GeneCaller.ncrnaFamilies.get(0).boundaryRefinerFactory instanceof CobalaminCmEndpointRefiner.Resources,
					"Main covariance selection must precede opening a declared but absent old tip configuration");
				final NcrnaFamilyLoader.Options conflict=new NcrnaFamilyLoader.Options();
				conflict.resources.root=source.getPath();
				conflict.boundaryOverride=true;
				boolean rejected=false;
				try{NcrnaFamilyLoader.load(FAMILIES[0], conflict, new FamilyResourceResolver());}
				catch(IllegalArgumentException expected){rejected=true;}
				check(rejected, "An explicit old main-tip request cannot silently override the packaged covariance choice");
			}finally{
				GeneCaller.ncrnaFamilies.clear();
				GeneCaller.ncrnaFamilies.addAll(previous);
				GeneCaller.cobalaminCallerResources=old;
				CallGenes.COBALAMIN_ENABLED=enabled;
				CallGenes.NCRNA_FAMILIES_ENABLED=ncrna;
				CallGenes.COBALAMIN_RESOURCE_ROOT=oldRoot;
			}
			pp.outstream.println("COBALAMIN_ORDINARY_PACKAGE_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void load(final File root){
		GeneCaller.ncrnaFamilies.clear();
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.COBALAMIN_ENABLED=true;
		CallGenes.COBALAMIN_RESOURCE_ROOT=root.getPath();
		CallGenes.loadCobalaminResources();
	}
	private static GeneCaller caller(){return new GeneCaller(60, 10, 10, -1, -1, -1, -1, -1, new GeneModel(false));}

	/** Constructor resets use real parsing with automatic family loading disabled.
	 * Explicit load(source) below must reach the conflict being tested, even when
	 * the preceding negative fixture left an incomplete root in static state. */
	private static void resetOptions(final File source, final File root){
		final File fasta=new File(root, "constructor.fa");
		try(Output out=new Output(fasta.getPath())){out.writer.print(">constructor\nACGTACGTACGTACGTACGT\n");}
		for(final String option:new String[]{"cobalamincm=off", "cobalamincmga=f", "cobalamincmmembers=all"}){
			try(CallGenes explicit=new CallGenes(new String[]{"in="+fasta.getPath(), "out=null", "taxonomy=f", "rna=f", "ncrna=f", "cobalamin=f", option})){
				boolean rejected=false;
				try{load(source);}catch(IllegalArgumentException expected){
					rejected=expected.getMessage()!=null && expected.getMessage().contains("The packaged cobalamin caller requires learned main-member covariance ends with GA hold");
				}
				check(rejected, "An explicitly incompatible current option must fail: "+option);
			}
			try(CallGenes omitted=new CallGenes(new String[]{"in="+fasta.getPath(), "out=null", "taxonomy=f", "rna=f", "ncrna=f", "cobalamin=f"})){
				load(source);
				check(GeneCaller.cobalaminCallerResources!=null,
					"A later constructor with omitted CM options must not inherit an earlier explicit-conflict flag: "+option);
			}
		}
	}

	/** The independently spelled development setup uses the frozen selected controls. */
	private static NcrnaClanScavenger control(final File root, final NcrnaFamily[] families, final NcrnaPathScoreConfig.Cell scores){
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides o=new FamilyResourceResolver.Overrides();
		o.root=root.getPath();
		final String[] paths=new String[3];
		final String[][] names=new String[3][];
		final byte[][][] bases=new byte[3][][];
		final NcrnaOverlapWindowBuilder.Caps[] caps=new NcrnaOverlapWindowBuilder.Caps[3];
		for(int member=0; member<3; member++){
			paths[member]=resolver.resolve(FAMILIES[member], "ddl_profiles", o).path;
			names[member]=families[member].modelNames;
			bases[member]=families[member].library;
			caps[member]=NcrnaOverlapWindowBuilder.Caps.load(resolver.resolve(FAMILIES[member], "key_bounds", o).path,
				MEMBERS[member], families[member].kmerSet.size(), 99);
			caps[member].validateKeys(families[member].kmerSet.toArray());
		}
		final NcrnaDdlProfileResource.ConsensusPair profiles=NcrnaDdlProfileResource.loadClanWithConsensusReplacement(paths, names, bases);
		final GeneCaller caller=caller();
		caller.setCobalaminClan(true, new NcrnaDdlProfileResource[]{profiles.memberOnly}, false);
		final NcrnaClanScavenger clan=(NcrnaClanScavenger)caller.ncrnaScavenger(0);
		clan.setModelPatience(4);
		clan.setWinnerRule(NcrnaMatchedBaseWinner.MODEL_SCORE, .58);
		clan.setScorePatience(100);
		clan.setConsensusPolicy(profiles, false);
		clan.setDdlOnlyAdmission(10);
		clan.setPathScoreCell(scores);
		clan.setMemberProposalPads(new int[]{380, 179, 123});
		clan.setOverlapWindows(true, caps);
		clan.setWindowCeilings(new int[]{467, 233, 199});
		clan.setCandidateAcceptance(new CobalaminCmAcceptance(new CobalaminCmAcceptance.Resources(families, root.getPath(),
			"filter_aligned", .60f, resolver.resolve(FAMILIES[0], "acceptance_policy", o).path, null, profiles.memberOnly, null, null, 25), clan, 0));
		return clan;
	}

	/** Independently bind the declared score artifact and check relative changes
	 * against the member bases. The saved export contains absolute coefficients. */
	private static NcrnaPathScoreConfig.Cell bindScores(final File root, final NcrnaFamily[] families, final String name){
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.root=root.getPath();
		final NcrnaPathScoreConfig table=new NcrnaPathScoreConfig(new FamilyResourceResolver().resolve(FAMILIES[0], "path_policy", overrides).path);
		table.bind(families);
		final NcrnaPathScoreConfig.Cell scores=table.cell(name);
		int changed=0;
		for(int member=0; member<families.length; member++){
			final float base=member==0 ? 2f : 1f;
			for(int model=0; model<families[member].modelNames.length; model++){
				check(scores.offset(member, model)==0, "Neither released comparison cell introduces an identity offset");
				if(scores.multiplier(member, model)!=base){changed++;}
			}
		}
		check(changed==(name.equals("greedy") ? 10 : 0),
			"The actual bound package must contain exactly the released zero or ten relative factor changes: "+name+" changed="+changed);
		return scores;
	}

	/** Select a real discovery-positive consensus; the greedy fixture requires
	 * a changed main-member coefficient so omitted factors affect score parity. */
	private static byte[] probe(final NcrnaFamily family, final NcrnaPathScoreConfig.Cell scores, final boolean changed){
		final long mask=(1L<<34)-1;
		for(int index=0; index<family.library.length; index++){
			if(changed && scores.multiplier(0, index)==2f){continue;}
			final byte[] model=family.library[index];
			if(model.length<150 || model.length>400){continue;}
			long word=0, first=-1;
			int valid=0;
			for(final byte base:model){
				final int n=base>=0 && base<AminoAcid.baseToNumber.length ? AminoAcid.baseToNumber[base] : -1;
				if(n<0){valid=0;word=0;continue;}
				word=((word<<2)|n)&mask;
				if(++valid<17 || !family.kmerSet.contains(word)){continue;}
				if(first<0){first=word;}
				else if(first!=word){
					final byte[] result=new byte[model.length+32];
					java.util.Arrays.fill(result, (byte)'T');
					System.arraycopy(model, 0, result, 16, model.length);
					return result;
				}
			}
		}
		throw new AssertionError("The serving library needs one bounded discovery-positive model, with a changed coefficient for the greedy fixture");
	}

	private static File cloneGroup(final File source, final File target, final String fault) throws Exception{
		check(target.mkdir(), "Each negative package owns a fresh data-only directory");
		final FamilyResourceResolver.Overrides o=new FamilyResourceResolver.Overrides();
		o.root=source.getPath();
		final String wrongProfile=new FamilyResourceResolver().resolve(FAMILIES[1], "ddl_profiles", o).path;
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(new File(source, "manifest.tsv").getPath()); Output out=new Output(new File(target, "manifest.tsv").getPath())){
			out.writer.println(input.reader.nextLine());
			out.writer.println(input.reader.nextLine());
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				final String member=p.parseString(0), role=p.parseString(1), path=p.parseString(4);
				if(fault.equals("missing_policy") && role.equals("caller_policy")){continue;}
				if(fault.equals("missing_profiles") && member.equals("RF00174") && role.equals("ddl_profiles")){continue;}
				out.writer.println(row);
				final String relative=member.equals(".") ? path : member+'/'+path;
				final File origin=fault.equals("wrong_profiles") && member.equals("RF00174") && role.equals("ddl_profiles")
					? new File(wrongProfile) : new File(source, relative);
				link(origin, new File(target, relative));
				if(role.equals("positional_config")){
					try(ConfigInput config=new ConfigInput(new File(source, relative).getPath())){
						while(config.next()){
							if(payload(config.key)){link(new File(source, member+'/'+config.value), new File(target, member+'/'+config.value));}
						}
					}
				}
			}
			if(fault.equals("unused_main_tips")){out.writer.println("RF00174\tpositional_config\t.\t.\tDO_NOT_OPEN.config\tconfig\t-\t-");}
		}
		return target;
	}
	private static void link(final File source, final File target) throws Exception{
		check(source.isFile(), "The synthetic package must link an existing declared payload: "+source);
		if(Files.exists(target.toPath())){return;}
		Files.createDirectories(target.toPath().getParent());
		Files.createSymbolicLink(target.toPath(), source.getAbsoluteFile().toPath());
	}
	private static boolean payload(final String key){
		return key.equals("library") || key.equals("traininghbm") || key.equals("modeloffsets")
			|| key.equals("net5") || key.equals("net3") || key.startsWith("table5") || key.startsWith("table3");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
}
