package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/** Root's finding, 2026-09-09: CellNetParser.parseHeader() defaults `dense=true` when neither
 * #dense nor #sparse appears in the header (ml/CellNetParser.java:334 field default, never
 * overwritten if no marker line matches) -- it does NOT throw. The prior parseDenseFlag() in
 * MagQCNetBundle.java threw "missing #dense/#sparse" instead, STRICTER than the real
 * producer/consumer contract: a legitimate net the production parser loads fine (defaulting to
 * dense) would have been wrongly REJECTED by pack()/validate(). This test constructs exactly that
 * case: takes a real legacy fixture net (which normally HAS an explicit #dense line), strips that
 * one line to simulate a net that never had it, packs a bundle containing the marker-less net
 * alongside the other 175 unmodified nets, and confirms (a) pack() now SUCCEEDS instead of
 * throwing, and (b) the marker-less net's score is BIT-IDENTICAL to the original (with-marker)
 * net's score on the same input -- proving the default correctly reconstructs dense=true, not
 * just that no exception was thrown. */
public final class MagQCNetBundleDefaultDenseTest {
	public static void main(String[] args) throws Exception {
		String aggManifest=null, subnetManifest=null, familyList=null, taxpgm=null,
			netroot=null, subsetroot=null, referenceBundlePath=null;
		for(String a:args) {
			int e=a.indexOf('='); if(e<1) throw new IllegalArgumentException("Expected key=value: "+a);
			String k=a.substring(0,e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("aggmanifest")) aggManifest=v;
			else if(k.equals("subnetmanifest")) subnetManifest=v;
			else if(k.equals("familylist")) familyList=v;
			else if(k.equals("taxpgm")) taxpgm=v;
			else if(k.equals("netroot")) netroot=v;
			else if(k.equals("subsetroot")) subsetroot=v;
			else if(k.equals("referencebundle")) referenceBundlePath=v;
		}
		if(aggManifest==null||subnetManifest==null||familyList==null||taxpgm==null||netroot==null||
				subsetroot==null||referenceBundlePath==null)
			throw new IllegalArgumentException("Required: aggmanifest= subnetmanifest= familylist= taxpgm= "+
				"netroot= subsetroot= referencebundle= (a legacy 176-subnet fixture's pack() inputs, plus "+
				"an already-packed bundle from the SAME unmodified fixture for score comparison; netroot "+
				"must contain a subnet id \"famset_0\" whose net carries an explicit #dense line)");

		Path origNet=Paths.get(netroot,"net_famset_0.bbnet");
		if(!Files.isRegularFile(origNet)) origNet=Paths.get(netroot,"nets","net_famset_0.bbnet");
		List<String> lines=Files.readAllLines(origNet, StandardCharsets.UTF_8);
		boolean hadDenseLine=false;
		List<String> stripped=new ArrayList<String>();
		for(String l:lines) {
			if(l.startsWith("#dense")) { hadDenseLine=true; continue; }
			stripped.add(l);
		}
		check(hadDenseLine,"test setup invalid: net_famset_0.bbnet has no #dense line to strip -- "+
			"fixture assumption violated, cannot exercise the default-when-absent path");
		for(String l:stripped) check(!l.startsWith("#sparse"),"test setup invalid: stripped net still has "+
			"an explicit #sparse line -- this would exercise the EXPLICIT-sparse path, not the "+
			"default-when-absent path");

		Path scratch=Files.createTempDirectory("magqc-default-dense-test-");
		Path modifiedNetRoot=scratch.resolve("nets_no_dense_marker");
		Files.createDirectories(modifiedNetRoot);
		Path realNetDir=Files.isDirectory(Paths.get(netroot,"nets")) ? Paths.get(netroot,"nets") : Paths.get(netroot);
		try(DirectoryStream<Path> ds=Files.newDirectoryStream(realNetDir)) {
			for(Path p:ds) Files.copy(p, modifiedNetRoot.resolve(p.getFileName()), StandardCopyOption.REPLACE_EXISTING);
		}
		Path modifiedNet=modifiedNetRoot.resolve("net_famset_0.bbnet");
		Files.write(modifiedNet, stripped, StandardCharsets.UTF_8);
		System.out.println("Stripped #dense line from net_famset_0.bbnet ("+lines.size()+" -> "+
			stripped.size()+" lines); rest of the 176-net fixture unchanged.");

		Map<String,String> packArgs=new LinkedHashMap<String,String>();
		packArgs.put("aggmanifest", aggManifest);
		packArgs.put("subnetmanifest", subnetManifest);
		packArgs.put("familylist", familyList);
		packArgs.put("taxpgm", taxpgm);
		packArgs.put("netroot", modifiedNetRoot.toString());
		packArgs.put("subsetroot", subsetroot);
		Path outBundle=scratch.resolve("marker_less_legacy.bbnets");
		packArgs.put("out", outBundle.toString());

		try {
			MagQCNetBundle.pack(packArgs);
		} catch(IOException e) {
			throw new AssertionError("pack() FAILED on a net with no explicit #dense/#sparse marker -- "+
				"the default-dense fix did not take effect: "+e, e);
		}
		System.out.println("pack() succeeded on a fixture containing a net with NO explicit #dense/#sparse "+
			"marker (previously this threw \"missing #dense/#sparse\").");

		MagQCNetBundle original=MagQCNetBundle.loadMultiOutput(Paths.get(referenceBundlePath));
		MagQCNetBundle markerLess=MagQCNetBundle.loadMultiOutput(outBundle);
		MagQCNetBundle.Subnet origSubnet=original.subnet("famset_0");
		MagQCNetBundle.Subnet markerLessSubnet=markerLess.subnet("famset_0");
		check(origSubnet.expectedInputs==markerLessSubnet.expectedInputs,"input width mismatch: orig="+
			origSubnet.expectedInputs+" markerLess="+markerLessSubnet.expectedInputs);

		float[] probe=new float[origSubnet.expectedInputs];
		for(int i=0;i<probe.length;i++) probe[i]=(float)((i*5)%11)/13f;
		float origScore=origSubnet.score(probe.clone());
		float markerLessScore=markerLessSubnet.score(probe.clone());
		check(Float.floatToIntBits(origScore)==Float.floatToIntBits(markerLessScore),
			"DEFAULT-DENSE FAILURE: marker-less net's score ("+markerLessScore+") does not bit-match the "+
			"original with-marker net's score ("+origScore+") on the same input -- the reconstructed dense "+
			"flag is wrong, not just 'no exception'");

		System.out.println("PASS: net with NO explicit #dense/#sparse marker packs successfully and scores "+
			"bit-identically ("+markerLessScore+") to the original with-marker net -- confirms the "+
			"reconstructed default is dense=TRUE, matching CellNetParser's real field default, not merely "+
			"'no exception thrown'.");
		System.out.println("MAGQC_BUNDLE_DEFAULT_DENSE PASS");
	}

	static void check(boolean ok,String what){ if(!ok) throw new AssertionError("FAIL: "+what); }
}
