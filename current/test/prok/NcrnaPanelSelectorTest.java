package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Fixed original-input panel selection; FASTA and GFF payloads do not exist.
 * Tests overlapping protected support, retained zeros, taxonomy balancing and
 * malformed/incomplete bindings without running a caller or changing reference.
 * @author Ganyu, Raiden
 */
public final class NcrnaPanelSelectorTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPanelSelectorTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){
					throw new IllegalArgumentException("Unknown parameter: "+arg);
				}
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append, "Supply out=fresh_fixture_directory only; no overwrite or input payload");
			new NcrnaPanelSelectorTest(new File(parser.out1).getAbsoluteFile()).run();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	private NcrnaPanelSelectorTest(final File root_){
		root=root_;
		check(root.mkdir(), "Use a fresh fixture directory under an existing parent");
	}

	private void run(){
		write("members.tsv", "family\tmember\treference_model\ntmRNA\tmajor\tRF00023\ntmRNA\tminor%A\tRF01849\ntmRNA\tminor%25B\tRF01850\n");
		write("policy.tsv", policy(5, 6, "controls"));
		final ByteBuilder census=new ByteBuilder("input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n");
		final ByteBuilder reference=new ByteBuilder("input_id\tfasta\tgff\n");
		final ByteBuilder memberships=new ByteBuilder("source\tmember_index\tinput_id\tsupplied_fasta\tsource_genome\tsource_assembly\ttaxid\tdeclared_compartment\treported_gff\treported_sha80\n");
		final ByteBuilder counts=new ByteBuilder(INPUT_HEADER);
		final ByteBuilder support=new ByteBuilder(MEMBER_HEADER);
		final ByteBuilder taxonomy=new ByteBuilder();
		for(int id=1; id<=8; id++){
			census.append(id).tab().append(fasta(id)).append("\t100\t1000\t1\n");
			reference.append(id).tab().append(fasta(id)).tab().append(path("unopened"+id+".gff")).nl();
			membership(memberships, "pool", id, id);
			final int scored=id<=6 ? 1 : 0;
			final int masked=id==8 ? 1 : 0;
			counts.append(id).append("\ttmRNA\t").append(scored).tab().append(masked).tab().append(scored>0)
				.tab().append(scored==0 && masked>0).tab().append(scored==0 && masked==0)
				.tab().append(id==1 ? 1 : 0).tab().append(id==2 ? 1 : 0).nl();
			for(int m=0; m<3; m++){
				final boolean positive=id<=6 && (id==2 ? m>0 : m==0);
				support.append(id).append("\ttmRNA\t");
				BuildNcrnaReference.field(support, NAMES[m]).tab();
				BuildNcrnaReference.field(support, MODELS[m]).tab().append(positive ? 1 : 0)
					.append("\t0\t").append(positive).append("\tNA\tNA\tNA\tNA\n");
			}
			if(id!=6){taxonomy.append("bacteria\t").append(id==2 || id==5 ? "B" : "A").append("\tg").append(id).nl();}
		}
		membership(memberships, "controls", 1, 7);
		membership(memberships, "controls", 2, 7); // Repeated source membership is still one control input.
		membership(memberships, "controls", 3, 8);
		write("census.tsv", census.toString());
		write("memberships.tsv", memberships.toString());
		write("reference.tsv", reference.toString());
		write("inventory.inputs.tsv", counts.toString());
		write("inventory.members.tsv", support.toString());
		provenance("inventory", "census.tsv");
		write("taxonomy.tsv", taxonomy.toString());
		write("taxonomy_manifest.tsv", taxonomyManifest("taxonomy.tsv"));
		select("selected", "inventory", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv");
		check(row("selected.summary.tsv", "tmRNA\t").startsWith("tmRNA\t8\t6\t5\t5\t1\t1\t7\t1\t"),
			"Five positive original inputs plus one masked-only and one zero control, including one unknown lineage");
		final ArrayList<String> config=lines("selected.inputs.config");
		final int[] expected={1,2,4,5,6,7,8};
		check(config.size()==expected.length, "Overlapping minority memberships and source aliases must not duplicate input IDs");
		for(int i=0; i<expected.length; i++){
			check(config.get(i).equals("in="+fasta(expected[i])), "Selected inputs retain canonical paths and numeric original-ID order");
		}
		check(row("selected.selection.tsv", "3\t").contains("NOT_SELECTED_BALANCED_LIMIT"),
			"TaxID10 already has mandatory input1, so phylumA filling prefers input4 with a new TaxID");
		check(row("selected.selection.tsv", "2\t").contains("PROTECTED_REFERENCE_MEMBER;MULTI_SUPPORTED"),
			"One original input preserves both minority member supports and multi-support");
		check(row("selected.support.tsv", "tmRNA\tminor%A\t").contains("\t1\t1\ttrue\t")
			&& row("selected.support.tsv", "tmRNA\tminor%25B\t").contains("\t1\t1\ttrue\t"),
			"Both escaped member names decode exactly once and retain all original support; a literal percent escape is not decoded twice");
		check(row("selected.strata.tsv", "tmRNA\tSOURCE\tcontrols\t").contains("\t2\t0\t2\t0\t"),
			"Repeated source membership is counted once in composition");
		check(lines("selected.reference.tsv").size()==8, "Every selected input retains its exact original reference pointer");
		check(!new File(fasta(1)).exists() && !new File(path("unopened1.gff")).exists(),
			"Metadata selection must not create or read missing genome/GFF payloads");

		// Input-table and taxonomy ordering are not an implicit RNG or tie-break.
		reverse("census.tsv", "reverse_census.tsv", true);
		reverse("inventory.inputs.tsv", "reverse.inputs.tsv", true);
		reverse("inventory.members.tsv", "reverse.members.tsv", true);
		provenance("reverse", "reverse_census.tsv");
		reverse("taxonomy.tsv", "reverse_taxonomy.tsv", false);
		write("reverse_taxonomy_manifest.tsv", taxonomyManifest("reverse_taxonomy.tsv"));
		select("reordered", "reverse", "reverse_census.tsv", "policy.tsv", "reverse_taxonomy_manifest.tsv");
		for(final String suffix:new String[]{"selection.tsv", "inputs.config", "reference.tsv", "summary.tsv", "support.tsv", "strata.tsv"}){
			check(lines("selected."+suffix).equals(lines("reordered."+suffix)), "Shuffling source rows must not change selected IDs or reports: "+suffix);
		}

		write("overflow_policy.tsv", policy(1, 1, "controls"));
		reject("overflow", "inventory", "census.tsv", "overflow_policy.tsv", "taxonomy_manifest.tsv", "Mandatory positive");
		write("large_policy.tsv", policy(7, 8, "controls"));
		reject("too_large", "inventory", "census.tsv", "large_policy.tsv", "taxonomy_manifest.tsv", "Insufficient positive");
		write("missing_source_policy.tsv", policy(5, 6, "absent"));
		reject("missing_source", "inventory", "census.tsv", "missing_source_policy.tsv", "taxonomy_manifest.tsv", "retained source");
		write("conflict_taxonomy.tsv", taxonomy.toString()+"bacteria\tDifferent\tg1\n");
		write("conflict_manifest.tsv", taxonomyManifest("conflict_taxonomy.tsv"));
		reject("conflict", "inventory", "census.tsv", "policy.tsv", "conflict_manifest.tsv", "Conflicting exact-source taxonomy");
		final ArrayList<String> badCounts=lines("inventory.inputs.tsv");
		badCounts.add(badCounts.get(1));
		writeLines("duplicate.inputs.tsv", badCounts);
		write("duplicate.members.tsv", support.toString());
		provenance("duplicate", "census.tsv");
		reject("duplicate_rejected", "duplicate", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv", "Repeated input/family");
		final ArrayList<String> missingSupport=lines("inventory.members.tsv");
		missingSupport.remove(1);
		write("missing.inputs.tsv", counts.toString());
		writeLines("missing.members.tsv", missingSupport);
		provenance("missing", "census.tsv");
		reject("missing_rejected", "missing", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv", "every declared inventory member");

		// These mutations keep each row's own positive flag consistent, but break
		// the cross-table invariant guaranteed by NcrnaSubtypeEvidence.support.
		final ArrayList<String> zeroUnion=lines("inventory.members.tsv");
		setSupport(zeroUnion, "7", "minor%25A", 1);
		write("zero_union.inputs.tsv", counts.toString());
		writeLines("zero_union.members.tsv", zeroUnion);
		provenance("zero_union", "census.tsv");
		reject("zero_union_rejected", "zero_union", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv", "exceeds the scored family union");
		final ArrayList<String> noSupport=lines("inventory.members.tsv");
		setSupport(noSupport, "3", "major", 0);
		write("no_support.inputs.tsv", counts.toString());
		writeLines("no_support.members.tsv", noSupport);
		provenance("no_support", "census.tsv");
		reject("no_support_rejected", "no_support", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv", "lacks sufficient qualified member support");
		new NcrnaPanelSelectorTest(new File(root, "exact_total")).exactTotal();
		System.err.println("NCRNA_PANEL_SELECTOR_TEST_PASS");
	}

	/** 1004 original IDs, overlapping member support and four unique controls.
	 * Original reference pointers retain both loci of input2 and the masked input.
	 */
	private void exactTotal(){
		write("members.tsv", "family\tmember\treference_model\ntmRNA\tmajor\tRF00023\ntmRNA\tminor%A\tRF01849\ntmRNA\tminor%25B\tRF01850\n");
		final ByteBuilder census=new ByteBuilder("input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n");
		final ByteBuilder reference=new ByteBuilder("input_id\tfasta\tgff\n");
		final ByteBuilder memberships=new ByteBuilder("source\tmember_index\tinput_id\tsupplied_fasta\tsource_genome\tsource_assembly\ttaxid\tdeclared_compartment\treported_gff\treported_sha80\n");
		final ByteBuilder counts=new ByteBuilder(INPUT_HEADER), support=new ByteBuilder(MEMBER_HEADER), taxonomy=new ByteBuilder();
		for(int id=1; id<=1004; id++){
			census.append(id).tab().append(fasta(id)).append("\t100\t1000\t1\n");
			reference.append(id).tab().append(fasta(id)).tab().append(path("unopened"+id+".gff")).nl();
			membership(memberships, "pool", id, id);
			final int scored=id>1001 ? 0 : id==2 ? 2 : 1;
			final int masked=id==1002 ? 1 : 0;
			counts.append(id).append("\ttmRNA\t").append(scored).tab().append(masked).tab().append(scored>0)
				.tab().append(scored==0 && masked>0).tab().append(scored==0 && masked==0).append("\t0\t").append(id==2 ? 1 : 0).nl();
			for(int m=0; m<3; m++){
				final boolean positive=id<=1001 && (m==0 ? id==1 || id>=201 : id==2 || (m==1 ? id>=3 && id<=101 : id>=102 && id<=200));
				support.append(id).append("\ttmRNA\t");
				BuildNcrnaReference.field(support, NAMES[m]).tab();
				BuildNcrnaReference.field(support, MODELS[m]).tab().append(positive ? scored : 0)
					.append("\t0\t").append(positive).append("\tNA\tNA\tNA\tNA\n");
			}
			taxonomy.append("bacteria\tP").append(id%3).append("\tg").append(id).nl();
		}
		membership(memberships, "controls", 1, 1);
		membership(memberships, "controls", 2, 2);
		membership(memberships, "controls", 3, 1002);
		membership(memberships, "controls", 4, 1003);
		membership(memberships, "controls", 5, 1003);
		write("census.tsv", census.toString());
		write("memberships.tsv", memberships.toString());
		write("reference.tsv", reference.toString());
		write("inventory.inputs.tsv", counts.toString());
		write("inventory.members.tsv", support.toString());
		provenance("inventory", "census.tsv");
		write("taxonomy.tsv", taxonomy.toString());
		write("taxonomy_manifest.tsv", taxonomyManifest("taxonomy.tsv"));
		write("policy.tsv", totalPolicy(1000, 100));
		select("selected", "inventory", "census.tsv", "policy.tsv", "taxonomy_manifest.tsv");
		final String summary=row("selected.summary.tsv", "tmRNA\t");
		check(summary.startsWith("tmRNA\t1004\t1001\t998\t998\t1\t1\t1000\t0\t")
			&& summary.endsWith("\t1000\t4\t2\t2\t996"),
			"Four original control IDs include two positives: 998 positive slots + two nonpositive controls = exactly1000, with996 additional positives");
		check(lines("selected.inputs.config").size()==1000 && lines("selected.reference.tsv").size()==1001,
			"Original-ID slot conservation must hold independently in caller and reference manifests");
		final String major=row("selected.support.tsv", "tmRNA\tmajor\t");
		final String minorA=row("selected.support.tsv", "tmRNA\tminor%A\t");
		final String minorB=row("selected.support.tsv", "tmRNA\tminor%25B\t");
		check(major.contains("\t802\t799\tfalse\t") && major.endsWith("\t300\t1\t802\t799"),
			"General fill retains799major-supported inputs after member quotas");
		check(minorA.contains("\t100\t100\tfalse\t") && minorA.endsWith("\t100\t2\t101\t101")
			&& minorB.contains("\t100\t100\tfalse\t") && minorB.endsWith("\t100\t3\t101\t101"),
			"Input2 credits both minority quotas once, retaining two scored loci per member without duplicate slots");
		check(row("selected.selection.tsv", "2\t").contains("RETAINED_SOURCE")
			&& row("selected.reference.tsv", "2\t").endsWith(path("unopened2.gff"))
			&& row("selected.reference.tsv", "1002\t").endsWith(path("unopened1002.gff")),
			"Supported fixed controls receive quota credit; full original reference and masked-control pointers remain unchanged");
		check(row("selected.selection.tsv", "1004\t").contains("NO_SCORED_TARGET_NOT_RETAINED_CONTROL"),
			"Nonpositive inputs outside the fixed controls cannot consume the positive fill budget");
		write("insufficient_policy.tsv", totalPolicy(1000, 101));
		reject("insufficient", "inventory", "census.tsv", "insufficient_policy.tsv", "taxonomy_manifest.tsv", "Insufficient supported inputs");
		write("budget_policy.tsv", totalPolicy(400, 100));
		reject("budget", "inventory", "census.tsv", "budget_policy.tsv", "taxonomy_manifest.tsv", "cannot fit the exact positive budget");
		write("unknown_policy.tsv", totalPolicy(1000, 100).replace("minor%A:100", "unknown:100"));
		reject("unknown", "inventory", "census.tsv", "unknown_policy.tsv", "taxonomy_manifest.tsv", "Unknown quota member");
		write("repeated_policy.tsv", totalPolicy(1000, 100).replace("minor%A:100", "major:100"));
		reject("repeated", "inventory", "census.tsv", "repeated_policy.tsv", "taxonomy_manifest.tsv", "Repeated quota member");
		reverse("census.tsv", "reverse_census.tsv", true);
		reverse("inventory.inputs.tsv", "reverse.inputs.tsv", true);
		reverse("inventory.members.tsv", "reverse.members.tsv", true);
		provenance("reverse", "reverse_census.tsv");
		select("reordered", "reverse", "reverse_census.tsv", "policy.tsv", "taxonomy_manifest.tsv");
		for(final String suffix:new String[]{"selection.tsv", "inputs.config", "reference.tsv", "summary.tsv", "support.tsv", "strata.tsv"}){
			check(lines("selected."+suffix).equals(lines("reordered."+suffix)), "Exact-total selection must be independent of input-table row ordering: "+suffix);
		}
	}

	private static String totalPolicy(final int total, final int minorityMinimum){
		return "family\tpositive_target\tpositive_ceiling\tprotected_members\tretained_sources\tkeep_protected_priors\tkeep_multi_supported\ttotal_inputs\tmember_min_inputs\n"
			+"tmRNA\t0\t"+total+"\t.\tcontrols\tfalse\tfalse\t"+total+"\tmajor:300,minor%A:"+minorityMinimum+",minor%25B:100\n";
	}

	/** Test-only corruption of one serialized row, retaining its local boolean. */
	private static void setSupport(final ArrayList<String> rows, final String input, final String encodedMember, final int scored){
		int changed=0;
		for(int i=1; i<rows.size(); i++){
			final String[] fields=rows.get(i).split("\t", -1);
			if(!fields[0].equals(input) || !fields[2].equals(encodedMember)){continue;}
			fields[4]=Integer.toString(scored);
			fields[6]=Boolean.toString(scored>0);
			rows.set(i, String.join("\t", fields));
			changed++;
		}
		check(changed==1, "Fixture mutation must alter exactly one original input/member row");
	}

	private void select(final String out, final String inventory, final String census, final String policy, final String taxonomy){
		BuildNcrnaReference.main(new String[]{"mode=selectpanel", "census="+path(census), "memberships="+path("memberships.tsv"),
			"reference="+path("reference.tsv"), "inventorymembers="+path("members.tsv"), "selectioninventory="+path(inventory),
			"selectionpolicy="+path(policy), "selectiontaxonomy="+path(taxonomy), "selectionfamily=tmRNA", "out="+path(out), "t=1"});
	}
	private void reject(final String out, final String inventory, final String census, final String policy, final String taxonomy, final String expected){
		try{
			select(out, inventory, census, policy, taxonomy);
			throw new AssertionError("Invalid metadata selection succeeded: "+out);
		}catch(IllegalArgumentException error){
			check(error.getMessage().contains(expected), "Wrong rejection: "+error);
		}
		check(!new File(path(out+".provenance.tsv")).exists(), "Rejected selection cannot publish completion");
		for(final String suffix:new String[]{".selection.tsv", ".inputs.config", ".reference.tsv"}){
			check(!new File(path(out+suffix)).exists(), "Rejected selection must fail before publishing a panel output: "+suffix);
		}
	}
	private void provenance(final String inventory, final String census){
		final ByteBuilder b=new ByteBuilder("role\tpath_or_value\n");
		for(final String[] pair:new String[][]{{"census", path(census)}, {"memberships", path("memberships.tsv")},
			{"reference", path("reference.tsv")}, {"members", path("members.tsv")},
			{"status", "REFERENCE_INVENTORY_COMPLETE_NO_SELECTION_OR_ADJUDICATION"}}){
			BuildNcrnaReference.field(b, pair[0]).tab();
			BuildNcrnaReference.field(b, pair[1]).nl();
		}
		write(inventory+".provenance.tsv", b.toString());
	}
	private String taxonomyManifest(final String map){
		return "schema\tpath\tsource_root\texpected_sha80\ngtdb_domain_phylum_basename_v1\t"+path(map)+"\t"+path("not_opened")+"\t00000000000000000000\n";
	}
	private static String policy(final int count, final int cap, final String source){
		return "family\tpositive_target\tpositive_ceiling\tprotected_members\tretained_sources\tkeep_protected_priors\tkeep_multi_supported\n"
			+"tmRNA\t"+count+"\t"+cap+"\tminor%A,minor%25B\t"+source+"\ttrue\ttrue\n";
	}
	private void membership(final ByteBuilder b, final String source, final int index, final int id){
		final String taxid=id==6 ? "NA" : id==1 || id==3 ? "10" : Integer.toString(id*10);
		b.append(source).tab().append(index).tab().append(id).tab().append(fasta(id)).append("\tg").append(id)
			.append("\tNA\t").append(taxid).append("\tprokaryotic\tNA\tNA\n");
	}
	private void reverse(final String original, final String out, final boolean header){
		final ArrayList<String> data=lines(original);
		final String first=header ? data.remove(0) : null;
		Collections.reverse(data);
		if(header){data.add(0, first);}
		writeLines(out, data);
	}
	private void writeLines(final String file, final ArrayList<String> data){
		final ByteBuilder b=new ByteBuilder();
		for(final String line:data){b.append(line).nl();}
		write(file, b.toString());
	}
	private void write(final String file, final String value){
		try(Output out=new Output(path(file))){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}
	}
	private ArrayList<String> lines(final String file){
		final ArrayList<String> values=new ArrayList<String>();
		try(Input in=new Input(path(file))){
			for(byte[] line; (line=in.reader.nextLine())!=null;){values.add(new String(line, StandardCharsets.UTF_8));}
		}
		return values;
	}
	private String row(final String file, final String prefix){
		String found=null;
		for(final String line:lines(file)){
			if(!line.startsWith(prefix)){continue;}
			check(found==null, "Repeated test result key: "+prefix);
			found=line;
		}
		check(found!=null, "Missing test result key: "+prefix);
		return found;
	}
	private String fasta(final int id){return path("not_opened/bacteria/g"+id+".fna.gz");}
	private String path(final String file){return new File(root, file).getAbsolutePath();}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private final File root;
	private static final String[] NAMES={"major", "minor%A", "minor%25B"}, MODELS={"RF00023", "RF01849", "RF01850"};
	private static final String INPUT_HEADER="input_id\tfamily\tscored_parents\tmasked_parents\tpositive_input\tmasked_only_input\tzero_target_input\tprotected_parents\tmulti_supported_scored_parents\n";
	private static final String MEMBER_HEADER="input_id\tfamily\tmember\treference_model\tqualified_scored_parents\tqualified_marginal_parents\tpositive_input\trecovery_scored_parents\trecovery_masked_parents\texclusive_scored_parents\texclusive_masked_parents\n";
}
