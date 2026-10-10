package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import parse.Parser;
import parse.LineParser1;
import fileIO.ByteFile;
import structures.ByteBuilder;

/** Compact-format fixtures and exact graph/metadata comparisons. @author Collei */
public final class HbmCompactTextProbe {

	public static void main(String[] args) throws Exception{
		String in=null, baseline=null, provenance=null, profile=null, reports=null;
		for(String arg:Parser.parseConfig(args)){
			if(arg.equalsIgnoreCase("selftest=t")){selfTest(); return;}
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected flag=value");}
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equalsIgnoreCase("in")){in=value;}
			else if(key.equalsIgnoreCase("baseline")){baseline=value;}
			else if(key.equalsIgnoreCase("provenance")){provenance=value;}
			else if(key.equalsIgnoreCase("profile")){profile=value;}
			else if(key.equalsIgnoreCase("reports")){reports=value;}
			else if(key.equalsIgnoreCase("t")){shared.Shared.setThreads(Integer.parseInt(value));}
			else if(key.equalsIgnoreCase("blockinput")){HbmCompactTextReader.BLOCK_INPUT=parse.Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown option "+key);}
		}
		if(reports!=null){
			final String[] paths=reports.split(",", -1);
			if(paths.length!=2){throw new IllegalArgumentException("reports=baseline.tsv,candidate.tsv required");}
			final Map<String,byte[]> first=reportRows(paths[0]), second=reportRows(paths[1]);
			if(!first.keySet().equals(second.keySet())){throw new AssertionError("Report row identities differ");}
			for(String key:first.keySet()){if(!Arrays.equals(first.get(key), second.get(key))){throw new AssertionError("Prediction differs: "+key);}}
			System.out.println("COMPACT_PROKCC_PREDICTIONS_PASS rows="+(first.size()-1)+" excluded_column=bin_worker_wall_seconds"); return;
		}
		if(in==null){throw new IllegalArgumentException("in=compact.hbmc required");}
		final HbmCompactTextReader.Result actual=HbmCompactTextReader.read(in);
		final ArrayList<String> ids=new ArrayList<String>(); final Map<String,byte[]> consensus=new HashMap<String,byte[]>();
		final HashMap<String,HbmCompactMetadata> expectedInfo=profile==null ? null : HbmCompactMetadata.fromProfile(profile);
		for(int i=0; i<actual.size(); i++){
			final String id=actual.name(i); ids.add(id); consensus.put(id, actual.consensus(i));
			if(expectedInfo!=null){equal(actual.metadata(i), expectedInfo.remove(id));}
		}
		if(expectedInfo!=null && !expectedInfo.isEmpty()){throw new AssertionError("Extra profile families");}
		if(baseline!=null){
			if(provenance==null){throw new IllegalArgumentException("provenance= required for native comparison");}
			final byte[][] pins=HbmBundleLoader.loadSemanticProvenance(provenance);
			for(int minCount:new int[]{1, 2, 10}){
				final HbmBundleLoader.Loaded expected=HbmBundleLoader.load(java.nio.file.Paths.get(baseline), ids, consensus::get, pins, minCount);
				final HbmBundleLoader.Loaded observed=HbmBundleLoader.load(java.nio.file.Paths.get(in), ids, consensus::get, pins, minCount);
				observed.assertStructuralEquivalent(expected);
			}
		}
		System.out.println("COMPACT_HBM_PASS families="+actual.size()+" filter="+actual.filter+
			" exact_graph_comparison="+(baseline!=null)+" exact_metadata_comparison="+(profile!=null));
	}

	/** Matches bins by path and compares all data fields except the declared timing column. */
	private static Map<String,byte[]> reportRows(String path){
		final Map<String,byte[]> result=new HashMap<String,byte[]>();
		final ByteFile in=ByteFile.makeByteFile(path, false); final LineParser1 row=new LineParser1('\t');
		int columns=-1;
		try{
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				row.set(line);
				if(row.termEquals("#columns", 0)){
					if(columns!=-1 || !row.termEquals("bin_worker_wall_seconds", row.terms()-1)){throw new IllegalArgumentException("Unexpected report columns");}
					columns=row.terms()-1; result.put("#columns", line);
				}else if(line.length>0 && line[0]!='#'){
					if(columns<1 || row.terms()!=columns){throw new IllegalArgumentException("Report data width");}
					final String key=row.parseString(1); row.length(columns-1);
					if(result.put(key, Arrays.copyOf(line, row.a()-1))!=null){throw new IllegalArgumentException("Duplicate bin path");}
				}
			}
		}finally{if(in.close()){throw new IllegalStateException("Report I/O error");}}
		if(result.size()<2){throw new IllegalArgumentException("Empty report");}
		return result;
	}

	private static void equal(HbmCompactMetadata a, HbmCompactMetadata b){
		if(b==null || a.identity!=b.identity || a.blosum!=b.blosum || a.hbm!=b.hbm ||
				a.start!=b.start || a.stop!=b.stop || a.minlen!=b.minlen || a.maxlen!=b.maxlen){
			throw new AssertionError("Metadata changed or family missing");
		}
	}

	/** Tests omission/partial headers, implicit topology, preserved sums and corrupt files. */
	private static void selfTest() throws Exception{
		final Path root=Files.createTempDirectory("hbm-compact-test-");
		final HbmCompactMetadata empty=new HbmCompactMetadata(), present=new HbmCompactMetadata();
		final ByteBuilder blank=new ByteBuilder(); empty.append(blank);
		if(blank.length()!=0){throw new AssertionError("Placeholder metadata was emitted");}
		present.identity=42.5; present.blosum=-12; present.hbm=0.8;
		present.start=0; present.stop=2; present.minlen=2; present.maxlen=9;
		final String good=fixture(present, false);
		final Path file=root.resolve("model.hbmc"); seal(file, good);
		final HbmCompactTextReader.Result result=HbmCompactTextReader.read(file.toString());
		if(result.size()!=1 || !Arrays.equals(result.consensus(0), Blosum62.encode("ARN".getBytes(StandardCharsets.US_ASCII), "fixture"))){
			throw new AssertionError("Consensus/roster changed");
		}
		equal(result.metadata(0), present);
		final AAGraph expected=new AAGraph(result.consensus(0), 0);
		for(int i=0; i<3; i++){
			Arrays.fill(expected.ref[i].count, 0); Arrays.fill(expected.ref[i].weight, 0);
			expected.ref[i].countSum=expected.ref[i].weightSum=10;
			expected.ref[i].count[result.consensus(0)[i]]=expected.ref[i].weight[result.consensus(0)[i]]=10;
		}
		expected.del[1].countSum=expected.del[1].weightSum=2;
		final AAGraphNode ins=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, 1);
		ins.countSum=ins.weightSum=2; ins.count[0]=ins.weight[0]=2; expected.ref[0].insEdge=ins;
		result.models().assertStructuralMatch(0, expected);
		final HbmCompactMetadata partial=new HbmCompactMetadata(); partial.hbm=0.25;
		final Path partialFile=root.resolve("partial.hbmc"); seal(partialFile, fixture(partial, true));
		equal(HbmCompactTextReader.read(partialFile.toString()).metadata(0), partial);
		final String[] changes={
			good.replace("#alphabet\tARN", "#alphabet\tAAN"),
			good.replace("#rows\t4", "#rows\t5"),
			good.replace("#start\t0", "#start\t-1"),
			good.replace("#stop\t2", "#stop\t3"),
			good.replace("#minlen\t2", "#minlen\t10"),
			good.replace("#minlen\t2", "#minlen\t2\n#minlen\t2"),
			good.replace("hbm=0.8", "hbm=NaN"),
			good.replace("hbm=0.8", "hbm=-999999"),
			good.replace("+\t2\t2", "+\t2\t3"),
			good.replace("+\t2\t2", "+\t2\t02"),
			good.replace("#filter\tnone", "#filter\tunknown"),
			good.replace("#rows\t4", "#unknown\t1\n#rows\t4")
		};
		for(int i=0; i<changes.length; i++){
			if(good.equals(changes[i])){throw new AssertionError("Unchanged fixture mutation");}
			final Path bad=root.resolve("bad"+i+".hbmc"); seal(bad, changes[i]); rejected(bad);
		}
		final byte[] bytes=Files.readAllBytes(file);
		final Path corrupt=root.resolve("checksum.hbmc");
		bytes[bytes.length-3]=(byte)(bytes[bytes.length-3]=='0' ? '1' : '0'); Files.write(corrupt, bytes); rejected(corrupt);
		final Path truncated=root.resolve("truncated.hbmc"); Files.write(truncated, good.getBytes(StandardCharsets.UTF_8)); rejected(truncated);
		// Standalone reading survives ordinary gzip as well as raw input.
		final Path zipped=root.resolve("model.hbmc.gz");
		try(java.util.zip.GZIPOutputStream out=new java.util.zip.GZIPOutputStream(Files.newOutputStream(zipped))){out.write(Files.readAllBytes(file));}
		HbmCompactTextReader.read(zipped.toString()).models().assertStructuralEquivalent(result.models());
		final Path crlf=root.resolve("crlf.hbmc"), noLf=root.resolve("no-final-lf.hbmc"), changed=root.resolve("changed.hbmc");
		final String sealed=new String(Files.readAllBytes(file), StandardCharsets.UTF_8);
		Files.write(crlf, sealed.replace("\n", "\r\n").getBytes(StandardCharsets.UTF_8));
		Files.write(noLf, sealed.substring(0, sealed.length()-1).getBytes(StandardCharsets.UTF_8));
		HbmCompactTextReader.read(crlf.toString()).models().assertStructuralEquivalent(result.models());
		HbmCompactTextReader.read(noLf.toString()).models().assertStructuralEquivalent(result.models());
		Files.write(changed, sealed.replace(":\t:\n", ";\t;\n").getBytes(StandardCharsets.UTF_8)); rejected(changed);
		// A header and body each cross multiple64KiB buffer boundaries, under LF and CRLF.
		final int length=70000; final char[] bases=new char[length]; Arrays.fill(bases, 'A');
		final String prefix=good.substring(0, good.indexOf("#name\t"));
		final ByteBuilder big=new ByteBuilder(prefix);
		big.append("#name\tlarge\n#alphabet\tA\n#consensus\t").append(new String(bases)).append("\n#rows\t").append(length).nl();
		for(int i=0; i<length; i++){big.append(":\t:\n");}
		final Path large=root.resolve("large.hbmc"), largeCrlf=root.resolve("large-crlf.hbmc");
		seal(large, big.toString());
		Files.write(largeCrlf, new String(Files.readAllBytes(large), StandardCharsets.UTF_8).replace("\n", "\r\n").getBytes(StandardCharsets.UTF_8));
		HbmCompactTextReader.read(large.toString()).models().assertStructuralEquivalent(HbmCompactTextReader.read(largeCrlf.toString()).models());
		System.out.println("COMPACT_HBM_FIXTURE_PASS omitted/partial headers, topology, lossy sums, gzip, malformed input, checksum");
		// These tiny disposable fixtures have no external dependencies.
		for(int i=0; i<changes.length; i++){Files.delete(root.resolve("bad"+i+".hbmc"));}
		for(Path path:Arrays.asList(file, partialFile, corrupt, truncated, zipped, crlf, noLf, changed, large, largeCrlf)){Files.delete(path);}
		Files.delete(root);
	}

	private static String fixture(HbmCompactMetadata info, boolean lossy){
		final ByteBuilder text=new ByteBuilder();
		text.append("#format\t").append(HbmCompactTextReader.FORMAT).nl();
		text.append("#graph_contract\t").append(HbmDenseTextBundle.CONTRACT).nl();
		text.append("#cutoff_units\t").append(HbmCompactTextReader.CUTOFF_UNITS).nl();
		text.append("#coordinates\tzero_based_inclusive\n#sums\toriginal\n#filter\t").append(lossy ? "singletons" : "none").nl();
		text.append("#columns\t+\tsum\tcounts\t-deletion\n#n_families\t1\n");
		for(int i=0; i<16; i++){text.append("#provenance_").append(i).append("\t11111111111111111111\n");}
		text.append("#name\tfixture\n#alphabet\tARN\n#consensus\tARN\n"); info.append(text);
		text.append("#rows\t4\n:\t").append(lossy ? "9" : ":").append("\n+\t2\t2\n:\t0\t:\t-2\n:\t0\t0\t:\n");
		return text.toString();
	}
	private static void seal(Path file, String text){
		final HbmCompactTextReader.Output out=new HbmCompactTextReader.Output(file.toString());
		try{out.print(text); out.finish();}finally{if(out.close()){throw new IllegalStateException("Fixture write failure");}}
	}
	private static void rejected(Path path) throws Exception{
		try{HbmCompactTextReader.read(path.toString());}
		catch(IllegalArgumentException | IOException e){return;}
		throw new AssertionError("Malformed compact file was accepted: "+path);
	}
	private HbmCompactTextProbe(){}
}
