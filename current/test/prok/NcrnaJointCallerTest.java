package test.prok;

import static test.prok.NcrnaCandidateCompetitionTest.equalGff;
import static test.prok.NcrnaCandidateCompetitionTest.metrics;
import static test.prok.NcrnaCandidateCompetitionTest.run;
import static test.prok.NcrnaCandidateCompetitionTest.write;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.FamilyResourceResolver;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaJointPolicy;
import prok.NcrnaProcess;
import prok.NcrnaReferenceCalls;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaPolicySearch.Measurement;
import prok.NcrnaSearchCoverage;
import prok.NcrnaStageDiagLogSink;
import prok.NcrnaTextIO.Input;
import prok.RrnaResourceIO;
import shared.Shared;

/** Bounded real-caller equivalence for both complete grouped families.
 * Ordinary controls are generated afresh, per original input. All policy GFFs
 * are independently regraded with the native saved-call grader. No optimization
 * is performed. The input recipe supplies qualified positive controls and pins.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaJointCallerTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaJointCallerTest.class, false);
		try{
			final Parser parser=new Parser();
			String reference=null, sixsMap=null, cobalaminMap=null, membersPath=null, lociPath=null;
			String policyControls=null, sam1Map=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("reference")){reference=b;}
				else if(a.equals("sixsmap")){sixsMap=b;}
				else if(a.equals("cobalaminmap")){cobalaminMap=b;}
				else if(a.equals("members")){membersPath=b;}
				else if(a.equals("loci")){lociPath=b;}
				else if(a.equals("policycontrols")){policyControls=b;}
				else if(a.equals("sam1map")){sam1Map=b;}
				else{check(parser.parse(arg, a, b), "Unknown test parameter: "+arg);}
			}
			check(parser.in1!=null && parser.out1!=null && reference!=null && sixsMap!=null && cobalaminMap!=null,
				"Supply in=cases.tsv reference=registered.tsv sixsmap=rules cobalaminmap=rules out=fresh_directory");
			check((membersPath==null)==(lociPath==null), "Optional native member coverage needs both members= and loci=");
			check(policyControls==null || sam1Map!=null, "Ordinary-policy controls require the SAM-I grading map");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(!root.exists() && root.mkdir(), "Fixture needs a fresh directory under an existing parent");
			final ArrayList<Case> cases=load(parser.in1, reference, policyControls==null ? 4 : 6);
			if(policyControls!=null){
				ordinaryPolicies(root, "6S", new String[]{"sixs_rf00013", "sixs_rf01685"}, sixsMap, reference, cases, policyControls);
				ordinaryPolicies(root, "cobalamin", new String[]{"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"}, cobalaminMap, reference, cases, policyControls);
				ordinaryPolicies(root, "sam1", new String[]{"sam1"}, sam1Map, reference, cases, policyControls);
				write(new File(root, "result.txt"), new String[]{"ORDINARY_POLICY_TABLE_CALLS_AND_LOGICAL_WORK_MATCH_SWEEP"});
				return;
			}
			group(root, "6S", new String[]{"sixs_rf00013", "sixs_rf01685"}, 60, sixsMap, reference, cases, membersPath, lociPath);
			group(root, "cobalamin", new String[]{"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"}, 100, cobalaminMap, reference, cases, membersPath, lociPath);
			write(new File(root, "result.txt"), new String[]{"JOINT_CALLER_GFF_GRADE_AND_LOGICAL_WORK_EQUIVALENCE_PASS"});
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Tests the public ordinary CLI against native sweep tuples, without a search.
	 * Each group gets exactly two original inputs and baseline/selected policies.
	 * Complete GFF equality includes every competing caller, not only the target.
	 */
	private static void ordinaryPolicies(final File root, final String family, final String[] members,
			final String map, final String reference, final ArrayList<Case> all, final String policyDirectory){
		final File directory=new File(root, family);
		check(directory.mkdir(), "Ordinary-policy fixture directory must be fresh");
		final ArrayList<Case> cases=new ArrayList<Case>();
		for(final Case c:all){if(c.family.equals(family)){cases.add(c);}}
		check(cases.size()==2, "Each ordinary-policy family needs two original validation inputs");
		final File policyFile=new File(policyDirectory, family+".tsv");
		final NcrnaJointPolicy[] points=NcrnaJointPolicy.load(policyFile.getPath());
		check(points.length==2, "Supply exactly shipped and selected tuples for each ordinary-policy group");
		final long[] ordinaryWork=new long[2];
		final ArrayList<String> common=new ArrayList<String>(Arrays.asList(COMMON));
		for(final Case c:cases){common.add("in="+c.fasta);}
		write(new File(directory, "caller.config"), common.toArray(new String[common.size()]));
		for(int p=0; p<points.length; p++){
			check(points[p].size()==members.length, "Explicit policy must cover every expected target member");
			final ArrayList<String> rows=new ArrayList<String>();
			rows.add(NcrnaJointPolicy.HEADER);
			for(int m=0; m<members.length; m++){
				final NcrnaJointPolicy.Member member=points[p].member(m);
				check(member.name.equals(members[m]), "Policy tuple order must match the registered target members");
				rows.add(points[p].name+"\t"+member.name+"\t"+(member.adaptive ? "adaptive" : "fixed")
					+"\t"+member.topN+"\t"+member.floor+"\t"+member.minHits+"\t"+member.topFrac+"\t"+member.qFrac+"\t"+member.patience);
			}
			final File singlePolicy=new File(directory, "ordinary_"+p+".policy.tsv");
			write(singlePolicy, rows.toArray(new String[rows.size()]));
			for(final Case c:cases){
				final String stem="ordinary_"+c.id+"_"+p;
				final File config=new File(directory, stem+".config"), profile=new File(directory, stem+".cpu.tsv");
				write(config, COMMON, "in="+c.fasta, "out="+new File(directory, stem+".gff"), "ncrnapolicy="+singlePolicy);
				profiled("callgenes.sh", config, profile);
				for(final String member:members){ordinaryWork[p]+=column(new File(profile+".kernels.tsv"), member, "logical_alignments");}
				if(p==0){
					final File omitted=new File(directory, "omitted_"+c.id+".config");
					write(omitted, COMMON, "in="+c.fasta, "out="+new File(directory, "omitted_"+c.id+".gff"));
					profiled("callgenes.sh", omitted, new File(directory, "omitted_"+c.id+".cpu.tsv"));
					equalGff(new File(directory, "omitted_"+c.id+".gff"), new File(directory, stem+".gff"));
				}
			}
		}
		final File prefix=new File(directory, "sweep"), config=new File(directory, "sweep.config");
		write(config, new String[]{"caller="+new File(directory, "caller.config"), "reference="+reference,
			"familymap="+map, "target="+family, "grading="+(members.length>1 ? "union" : "target"),
			"joint=t", "policies="+policyFile, "out="+prefix, "savegff=t", "cachewindows=0", "reusecandidates=f"});
		profiled("dev/sweepncrnashortlist.sh", config, new File(directory, "sweep.cpu.tsv"));
		for(int p=0; p<points.length; p++){
			for(int i=0; i<cases.size(); i++){
				if(p==0){
					final Case c=cases.get(i);
					check(metrics(new File(prefix+".inputs.tsv"), c.id, points[p].name)[0]==c.expectedTp,
						"The shipped policy must retain each declared positive control: "+family+" / "+c.id);
				}
				equalGff(new File(directory, "ordinary_"+cases.get(i).id+"_"+p+".gff"),
					new File(prefix+".input"+i+".policy"+p+".gff"));
			}
			check(ordinaryWork[p]>0 && ordinaryWork[p]==column(new File(prefix+".summary.tsv"), points[p].name, "kernel_logical_alignments"),
				"Ordinary target-member logical work must equal the native sweep sum for "+family+" / "+points[p].name);
		}
		check(column(new File(prefix+".summary.tsv"), points[0].name, "TP")>0,
			"The ordinary-policy control must retain positive reference detections; empty equality is insufficient");
		write(new File(directory, "work.tsv"), new String[]{"policy\tordinary_logical_alignments",
			points[0].name+"\t"+ordinaryWork[0], points[1].name+"\t"+ordinaryWork[1]});
	}

	private static void group(final File root, final String family, final String[] members, final int cap,
			final String map, final String reference, final ArrayList<Case> all, final String membersPath, final String lociPath){
		final File directory=new File(root, family);
		check(directory.mkdir(), "Group fixture directory must be fresh");
		final ArrayList<Case> cases=new ArrayList<Case>();
		for(final Case c:all){if(family.equals(c.family)){cases.add(c);}}
		check(cases.size()==2, "Each bounded group fixture requires exactly two original inputs");
		final Rule[] rules=GradeNcrnaReference.loadRules(map);
		final long[] memberCalls=new long[rules.length], memberWork=new long[members.length];
		final ArrayList<String> common=new ArrayList<String>(Arrays.asList(COMMON));
		for(final Case c:cases){common.add("in="+c.fasta);}
		write(new File(directory, "caller.config"), common.toArray(new String[common.size()]));
		for(final Case c:cases){
			final File gff=new File(directory, "ordinary_"+c.id+".gff"), config=new File(directory, "ordinary_"+c.id+".config");
			write(config, COMMON, "in="+c.fasta, "out="+gff);
			profiled("callgenes.sh", config, new File(directory, "ordinary_"+c.id+".cpu.tsv"));
			final NcrnaReferenceCalls.Result calls=new NcrnaReferenceCalls(rules).read(gff.getPath(), c.id);
			for(int r=0; r<rules.length; r++){memberCalls[r]+=calls.ruleCounts[r];}
			for(int m=0; m<members.length; m++){
				memberWork[m]+=column(new File(directory, "ordinary_"+c.id+".cpu.tsv.kernels.tsv"), members[m], "logical_alignments");
			}
			final long[] grade=grade(directory, "ordinary_"+c.id, map, reference, new Case[]{c}, new File[]{gff}, family);
			check(grade[0]==c.expectedTp && grade[1]==0 && grade[2]==0,
				"Ordinary fixture must retain its declared qualified positive coverage: "+c.id+" "+Arrays.toString(grade));
		}
		for(int r=0; r<rules.length; r++){
			check(!rules[r].include || memberCalls[r]>0, "A requested member must actually emit calls; empty equality is insufficient: "+rules[r].modelPrefix);
		}
		for(int m=0; m<members.length; m++){check(memberWork[m]>0, "No ordinary verification work for intended member "+members[m]);}
		final ArrayList<String> policies=new ArrayList<String>();
		policies.add(NcrnaJointPolicy.HEADER);
		final NcrnaJointPolicy.Member[] packaged=new NcrnaJointPolicy.Member[members.length];
		for(int i=0; i<members.length; i++){packaged[i]=packagedMember(members[i]);}
		final String[] points={"baseline", "legacy_fixed", "legacy_adaptive", "changed", "first", "full"};
		for(final String point:points){
			for(final NcrnaJointPolicy.Member member:packaged){
				if(point.equals("baseline")){policies.add(policyRow(point, member));}
				else{
					final boolean changed=point.equals("changed"), fixed=point.equals("legacy_fixed");
					final int patience=point.equals("full") ? Integer.MAX_VALUE : point.equals("first") || changed ? 0 : -1;
					policies.add(policyRow(point, new NcrnaJointPolicy.Member(member.name, !fixed,
						changed ? 1 : cap, fixed ? 0f : changed ? 12f : 1f, 1,
						changed ? .75f : 0f, changed ? .03125f : 0f, patience)));
				}
			}
		}
		final int baseline=pointIndex(points, "baseline"), legacyFixed=pointIndex(points, "legacy_fixed");
		final int legacyAdaptive=pointIndex(points, "legacy_adaptive"), changed=pointIndex(points, "changed");
		final int first=pointIndex(points, "first"), full=pointIndex(points, "full");
		write(new File(directory, "policies.tsv"), policies.toArray(new String[policies.size()]));
		long[] work=null;
		for(int mode=0; mode<4; mode++){
			final String name="control"+mode;
			final File prefix=new File(directory, name), config=new File(directory, name+".config");
			write(config, new String[]{"caller="+new File(directory, "caller.config"), "reference="+reference,
				"familymap="+map, "target="+family, "grading=union", "joint=t", "policies="+new File(directory, "policies.tsv"),
				"out="+prefix, "savegff=t", "cachewindows="+((mode&1)==0 ? 0 : 128), "reusecandidates="+((mode&2)!=0)});
			profiled("dev/sweepncrnashortlist.sh", config, new File(directory, name+".cpu.tsv"));
			final long[] measured=new long[points.length];
			for(int p=0; p<points.length; p++){
				final String point=points[p];
				final File[] gffs=new File[cases.size()];
				for(int i=0; i<gffs.length; i++){
					gffs[i]=new File(prefix+".input"+i+".policy"+p+".gff");
					if(p==baseline){equalGff(new File(directory, "ordinary_"+cases.get(i).id+".gff"), gffs[i]);}
					if(p==legacyAdaptive || p==full){
						equalGff(new File(prefix+".input"+i+".policy"+legacyFixed+".gff"), gffs[i]);
					}
					if(mode>0){equalGff(new File(directory, "control0.input"+i+".policy"+p+".gff"), gffs[i]);}
				}
				final long[] actualGrade=grade(directory, name+"_"+point, map, reference,
					cases.toArray(new Case[cases.size()]), gffs, family);
				check(Arrays.equals(actualGrade, metrics(new File(prefix+".summary.tsv"), point)), "Saved-GFF and live native grades differ: "+name+" "+point);
				for(int i=0; i<gffs.length; i++){
					final Case c=cases.get(i);
					final long[] individual=grade(directory, name+"_"+point+"_input"+i, map, reference,
						new Case[]{c}, new File[]{gffs[i]}, family);
					check(Arrays.equals(individual, metrics(new File(prefix+".inputs.tsv"), c.id, point)),
						"Each original input's ten native grade fields must agree independently: "+c.id+" "+point);
				}
				measured[p]=column(new File(prefix+".summary.tsv"), point, "kernel_logical_alignments");
			}
			check(measured[legacyFixed]==measured[legacyAdaptive] && measured[full]==measured[legacyFixed],
				"Fixed1/adaptive1 and legacy/MAX full walks use identical mapping and must retain logical work");
			check(measured[changed]<measured[legacyFixed],
				"The declared positive-input witness must demonstrate reduced work for the restrictive legacy comparison");
			if(family.equals("cobalamin")){
				check(measured[first]<measured[legacyAdaptive],
					"The declared cobalamin positive witness must demonstrate less work at patience zero than the same full-walk mapping");
			}
			if(work==null){work=measured;}
			else{check(Arrays.equals(work, measured), "Cache/reuse modes must preserve summed logical work for every policy");}
		}
		long ordinaryWork=0;
		for(final long count:memberWork){ordinaryWork+=count;}
		check(work[baseline]==ordinaryWork, "Joint packaged baseline work must equal the sum of ordinary member kernel counts");
		stages(directory, family, map, reference, common, cases, points);
		final NcrnaJointPolicy.Member baselineMember=packaged[0];
		write(new File(directory, "legacy.tsv"), new String[]{"policy\ttopn\t"+(baselineMember.adaptive ? "floor" : "minhits"),
			"baseline\t"+baselineMember.topN+"\t"+(baselineMember.adaptive ? Float.toString(baselineMember.floor) : Integer.toString(baselineMember.minHits))});
		write(new File(directory, "legacy.config"), new String[]{"caller="+new File(directory, "caller.config"), "reference="+reference,
			"familymap="+map, "target="+members[0], "grading=union", "policies="+new File(directory, "legacy.tsv"),
			"out="+new File(directory, "legacy"), "savegff=t", "cachewindows=0", "reusecandidates=f"});
		profiled("dev/sweepncrnashortlist.sh", new File(directory, "legacy.config"), new File(directory, "legacy.cpu.tsv"));
		for(int i=0; i<cases.size(); i++){equalGff(new File(directory, "ordinary_"+cases.get(i).id+".gff"), new File(directory, "legacy.input"+i+".policy0.gff"));}
		if(membersPath!=null){coverage(directory, family, map, reference, common, cases, points, membersPath, lociPath);}
	}

	/** Two concurrent inputs have independent policy sinks and unchanged GFF/work. */
	private static void stages(final File directory, final String family, final String map, final String reference,
			final ArrayList<String> common, final ArrayList<Case> cases, final String[] points){
		final File caller=new File(directory, "stages.caller.config"), config=new File(directory, "stages.config");
		final ArrayList<String> args=new ArrayList<String>(common);
		args.add("t=2");
		write(caller, args.toArray(new String[args.size()]));
		final File prefix=new File(directory, "stages");
		write(config, new String[]{"caller="+caller, "reference="+reference, "familymap="+map,
			"target="+family, "grading=union", "joint=t", "policies="+new File(directory, "policies.tsv"),
			"out="+prefix, "savegff=t", "savestages=t", "cachewindows=128", "reusecandidates=t"});
		profiled("dev/sweepncrnashortlist.sh", config, new File(directory, "stages.cpu.tsv"));
		for(int i=0; i<cases.size(); i++){
			final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
			RrnaResourceIO.stream(cases.get(i).fasta, read -> {
				final String name=RrnaResourceIO.first(read.id);
				check(!lengths.contains(name), "Fixture input must have unique contig IDs for independent trace attribution");
				lengths.put(name, read.length());
			});
			for(int p=0; p<points.length; p++){
				final File old=new File(directory, "control3.input"+i+".policy"+p+".gff");
				equalGff(old, new File(prefix+".input"+i+".policy"+p+".gff"));
				check(!new File(directory, "control3.input"+i+".policy"+p+".stages.tsv").exists(),
					"Default-off tracing must create no stage sidecars");
				stageRows(new File(prefix+".input"+i+".policy"+p+".stages.tsv"), lengths);
			}
		}
		for(final String point:points){
			for(final String metric:new String[]{"TP", "FP", "FN", "kernel_logical_alignments", "aligned_model_bases"}){
				check(column(new File(prefix+".summary.tsv"), point, metric)==column(new File(directory, "control3.summary.tsv"), point, metric),
					"Enabling stage observation must preserve native detection and logical-work counts: "+metric);
			}
		}
		if(family.equals("cobalamin")){collisionStages(directory, map, points.length);}
	}

	/** Same header, distinct lengths/inputs: incorrect shared attribution must fail. */
	private static void collisionStages(final File directory, final String map, final int policies){
		final String[] names={"a", "b"};
		final ArrayList<String> caller=new ArrayList<String>(Arrays.asList(COMMON));
		final ArrayList<String> references=new ArrayList<String>();
		references.add("input_id\tfasta\tgff");
		for(int i=0; i<names.length; i++){
			final File fasta=new File(directory, "collision_"+names[i]+".fa"), gff=new File(directory, "collision_"+names[i]+".gff");
			final char[] bases=new char[500+i*100];
			Arrays.fill(bases, 'N');
			write(fasta, new String[]{">same_contig_header", new String(bases)});
			new NcrnaReferenceGff().write(gff.getPath(), names[i], new ArrayList<NcrnaReferenceGff.Row>());
			caller.add("in="+fasta);
			references.add(names[i]+"\t"+fasta+"\t"+gff);
		}
		caller.add("t=2");
		final File callerPath=new File(directory, "collision.caller.config"), referencePath=new File(directory, "collision.reference.tsv");
		write(callerPath, caller.toArray(new String[caller.size()]));
		write(referencePath, references.toArray(new String[references.size()]));
		for(int record=0; record<2; record++){
			final File prefix=new File(directory, "collision"+record), config=new File(directory, "collision"+record+".config");
			write(config, new String[]{"caller="+callerPath, "reference="+referencePath, "familymap="+map,
				"target=cobalamin", "grading=union", "joint=t", "policies="+new File(directory, "policies.tsv"),
				"out="+prefix, "savegff=t", "savestages="+(record==1), "cachewindows=128", "reusecandidates=t"});
			profiled("dev/sweepncrnashortlist.sh", config, new File(directory, "collision"+record+".cpu.tsv"));
		}
		for(int i=0; i<2; i++){
			final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
			lengths.put("same_contig_header", 500+i*100);
			for(int p=0; p<policies; p++){
				equalGff(new File(directory, "collision0.input"+i+".policy"+p+".gff"),
					new File(directory, "collision1.input"+i+".policy"+p+".gff"));
				stageRows(new File(directory, "collision1.input"+i+".policy"+p+".stages.tsv"), lengths);
			}
		}
	}

	/** Every emitted row belongs to the declared original input, even under t=2. */
	private static void stageRows(final File path, final ObjectIntMap<String> lengths){
		final LineParser1 parser=new LineParser1('\t');
		byte[] last=null;
		int length=-1;
		long rows=0;
		try(Input input=new Input(path.getPath())){
			final byte[] header=input.reader.nextLine();
			check(header!=null && new String(header, StandardCharsets.UTF_8).equals(NcrnaStageDiagLogSink.HEADER),
				"Bound stage files retain the existing diagnostic schema");
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				check(parser.terms()==21, "A stage row cannot mix buffers from two input/policy writers");
				parser.setBounds(2);
				boolean same=last!=null && last.length==parser.b()-parser.a();
				for(int j=0; same && j<last.length; j++){same=last[j]==line[parser.a()+j];}
				if(!same){
					last=parser.parseByteArray(2);
					length=lengths.get(RrnaResourceIO.first(new String(last, StandardCharsets.UTF_8)));
				}
				check(length>0 && length==parser.parseInt(4), "Every stage row must bind the exact input's contig and length");
				check(parser.termEquals("SEED", 0) || parser.termEquals("MODEL", 0) || parser.termEquals("WINDOW", 0)
					|| parser.termEquals("PATH", 0), "Unexpected stage event in a bound policy stream");
				rows++;
			}
		}
		check(rows>0, "The positive fixture must exercise the actual per-input/policy stage sink");
	}

	/** Opt-in accounting must preserve every complete GFF for each fixed tuple. */
	private static void coverage(final File directory, final String family, final String map, final String reference,
			final ArrayList<String> common, final ArrayList<Case> cases, final String[] points,
			final String membersPath, final String lociPath){
		final File prefix=new File(directory, "coverage"), config=new File(directory, "coverage.config");
		write(config, new String[]{"caller="+new File(directory, "caller.config"), "reference="+reference,
			"familymap="+map, "target="+family, "grading=union", "joint=t", "policies="+new File(directory, "policies.tsv"),
			"out="+prefix, "savegff=t", "cachewindows=0", "reusecandidates=f", "coverage=t",
			"members="+membersPath, "loci="+lociPath});
		profiled("dev/sweepncrnashortlist.sh", config, new File(directory, "coverage.cpu.tsv"));
		for(int p=0; p<points.length; p++){
			for(int i=0; i<cases.size(); i++){
				equalGff(new File(directory, "control0.input"+i+".policy"+p+".gff"), new File(prefix+".input"+i+".policy"+p+".gff"));
			}
		}
		final ArrayList<String> names=new ArrayList<String>();
		try(Input input=new Input(membersPath)){
			final LineParser1 parser=new LineParser1('\t');
			input.reader.nextLine();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.termEquals(family, 0)){names.add(parser.parseString(1));}
			}
		}
		check(names.size()==(family.equals("6S") ? 2 : 3), "Native coverage fixture requires the complete qualified reference-member roster");
		final Measurement[] measured=new Measurement[points.length];
		for(int p=0; p<points.length; p++){
			measured[p]=new Measurement(column(new File(prefix+".summary.tsv"), points[p], "TP"),
				column(new File(prefix+".summary.tsv"), points[p], "FP"),
				column(new File(prefix+".summary.tsv"), points[p], "kernel_logical_alignments"),
				column(new File(prefix+".summary.tsv"), points[p], "lost_baseline_TP"));
		}
		final NcrnaSearchCoverage.Counts[] counts=new NcrnaSearchCoverage(common.toArray(new String[common.size()]),
			reference, names.toArray(new String[names.size()])).read(prefix.getPath(), points, measured);
		final int fixed=pointIndex(points, "legacy_fixed"), adaptive=pointIndex(points, "legacy_adaptive");
		check(counts[fixed].lastLoss()==counts[adaptive].lastLoss() && counts[fixed].firstGain()==counts[adaptive].firstGain(),
			"Equivalent legacy mappings must have identical coverage changes relative to the packaged baseline");
		for(int m=0; m<names.size(); m++){
			check(counts[fixed].memberLastLoss(m)==counts[adaptive].memberLastLoss(m)
				&& counts[fixed].memberFirstGain(m)==counts[adaptive].memberFirstGain(m),
				"Equivalent legacy mappings must preserve each qualified member's coverage changes");
		}
	}

	/** Integration input follows the real package; independent literal defaults live in NcrnaCallerPolicyTest. */
	private static NcrnaJointPolicy.Member packagedMember(final String name){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.endpoints=false;
		options.boundaryOverride=Boolean.FALSE;
		final NcrnaFamily family=NcrnaFamilyLoader.load(name, options, new FamilyResourceResolver());
		return new NcrnaJointPolicy.Member(name, family.adaptive, family.indexTopN, family.adaptFloor,
			family.fixedMinHits, family.adaptTopFrac, family.adaptQFrac, family.modelPatience);
	}
	private static String policyRow(final String point, final NcrnaJointPolicy.Member member){
		return point+"\t"+member.name+"\t"+(member.adaptive ? "adaptive" : "fixed")+"\t"+member.topN+"\t"+member.floor
			+"\t"+member.minHits+"\t"+member.topFrac+"\t"+member.qFrac+"\t"+member.patience;
	}
	private static int pointIndex(final String[] points, final String name){
		for(int i=0; i<points.length; i++){if(name.equals(points[i])){return i;}}
		throw new IllegalArgumentException("Missing intended policy fixture: "+name);
	}

	private static long[] grade(final File directory, final String name, final String map, final String reference,
			final Case[] cases, final File[] gffs, final String family){
		final ArrayList<String> rows=new ArrayList<String>();
		rows.add("input_id\tfasta\tgff");
		for(int i=0; i<cases.length; i++){rows.add(cases[i].id+"\t"+cases[i].fasta+"\t"+gffs[i]);}
		write(new File(directory, name+".queries.tsv"), rows.toArray(new String[rows.size()]));
		final File config=new File(directory, name+".grade.config");
		write(config, new String[]{"in="+new File(directory, name+".queries.tsv"), "reference="+reference,
			"familymap="+map, "restrictreference=t", "out="+new File(directory, name+".grade")});
		run("dev/gradencrnareference.sh", config, 1024);
		return metrics(new File(directory, name+".grade.summary.tsv"), family);
	}

	private static void profiled(final String script, final File config, final File profile){
		final String path=Data.findPath("?"+script, false);
		check(path!=null, "Missing installed launcher "+script);
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("env",
			"JAVA_TOOL_OPTIONS=-XX:ActiveProcessorCount=4 -Dbbtools.diagnostic.ncrnaCpuOut="+profile,
			"bash", path, "config="+config.getPath()));
		NcrnaProcess.addJvmOptions(command, 4096);
		NcrnaProcess.run(command, null, "Bounded grouped-caller validation");
	}

	/** Metadata accounting only; no sequence payloads or biological inference here. */
	private static long column(final File path, final String key, final String column){
		final LineParser1 parser=new LineParser1('\t');
		int index=-1;
		long result=-1;
		try(Input input=new Input(path.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				if(index<0){
					for(int i=0; i<parser.terms(); i++){if(parser.termEquals(column, i)){index=i;}}
					check(index>=0, "Missing required counter "+column+" in "+path);
				}else if(parser.termEquals(key, 0)){
					check(result<0 && !parser.termEquals("NA", index), "Counter must be unique and instrumented");
					result=parser.parseLong(index);
				}
			}
		}
		check(result>=0, "Missing nonnegative counter for "+key+" in "+path);
		return result;
	}

	private static ArrayList<Case> load(final String path, final String reference, final int expectedCases){
		final ArrayList<GradeNcrnaReference.Binding> bindings=GradeNcrnaReference.loadManifest(reference);
		final ArrayList<Case> cases=new ArrayList<Case>();
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			check(header!=null, "Missing case header");
			parser.set(header);
			check(parser.terms()==4 && parser.termEquals("family", 0) && parser.termEquals("input_id", 1)
				&& parser.termEquals("fasta", 2) && parser.termEquals("expected_TP", 3), "Case header must declare family, original input ID, FASTA and expected TP");
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				check(parser.terms()==4, "Invalid bounded case row");
				final Case c=new Case(parser.parseString(0), parser.parseString(1), parser.parseString(2), parser.parseInt(3));
				int matches=0;
				for(final GradeNcrnaReference.Binding binding:bindings){if(binding.id.equals(c.id) && binding.fasta.equals(c.fasta)){matches++;}}
				check(matches==1 && new File(c.fasta).isFile(), "Case must bind one registered original input and readable FASTA: "+c.id);
				cases.add(c);
			}
		}
		check(cases.size()==expectedCases, "Bounded recipe requires exactly "+expectedCases+" original inputs");
		return cases;
	}
	private static final class Case{
		Case(final String family_, final String id_, final String fasta_, final int tp_){
			family=family_;
			id=id_;
			fasta=fasta_;
			expectedTp=tp_;
		}
		final String family, id, fasta;
		final int expectedTp;
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] COMMON={"taxonomy=f", "callcds=t", "calltrna=t", "call16s=t", "call23s=t", "call5s=t",
		"call18s=f", "ncrna=t", "tmrna=t", "sixs=t", "fmn=t", "sam1=t", "tpp=t", "tppcdsstop=t", "glycine=t",
		"cobalamin=t", "euk=f", "ordered=t", "t=1", "ow=f", "append=f"};
}
