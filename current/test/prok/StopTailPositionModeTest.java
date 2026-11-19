package test.prok;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilyProfiles;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import prok.StopTailCascadeResources;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyTailCascade;
import prot.ProteinSequence;
import prot.ProteinSpectrum;
import structures.ByteBuilder;

/** Tests the caller's position selector against real profiles and query workers.
 * Structural fixture identities avoid hashing newly written test payloads.
 * @author Sayu
 */
public final class StopTailPositionModeTest{

	public static void main(String[] args)throws Exception{
		if(args.length<3 || args.length>4 || args.length==4 && !args[3].equals("legacy")){
			throw new IllegalArgumentException("Expected fresh directory, retained scalar table, "
				+"inherited scalar identity, optional legacy");
		}
		final Path root=Paths.get(args[0]).toAbsolutePath(); Files.createDirectory(root);
		final Path bounds=Paths.get(args[1]).toRealPath();
		final String boundsIdentity=args[2];
		check(Files.isRegularFile(bounds) && boundsIdentity.matches("[0-9a-f]{20}"),
			"The retained scalar loader requires an actual resource with its inherited identity");
		final byte[] query=("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"+
			"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY").getBytes(StandardCharsets.US_ASCII);
		check(query.length==100,"Fixture stays within the released 100-aa limit");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.exact(5);
		final char[] plain=new char[1024]; DDLFamilyProfiles.sketch(query,geometry,0,plain,null);
		final ProteinSpectrum.Spectrum spectrum=new ProteinSpectrum.Spectrum(2); spectrum.add(query);
		final FamilyReference family=new FamilyReference("q",0,1,100,96,new ProteinSequence("q",query),plain,
			FamilyReferenceExport.frequencies(spectrum),0.5f,50,0.5f,1.1f,99,1);
		final ByteBuilder refs=new ByteBuilder();
		FamilyReferenceCodec.appendHeader(new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,
			5,1024,0,1,1,0,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING),refs);
		FamilyReferenceCodec.appendFiveLines(family,refs);
		final Path references=write(root.resolve("references.tsv"),refs.toString());
		final char[] fixed=sketch(query,geometry,ProteinPositionMode.FIXED);
		final char[] quadrant=sketch(query,geometry,ProteinPositionMode.FIXED_QUADRANT);
		check(!Arrays.equals(fixed,quadrant),"The two layouts must differ so fallback cannot pass this fixture");
		final Path fixedProfiles=profiles(root,"fixed",fixed,ProteinPositionMode.FIXED);
		final Path quadrantProfiles=profiles(root,"quadrant",quadrant,ProteinPositionMode.FIXED_QUADRANT);
		final String settings="familyreferences\t"+references+"\nfamilyreferencessha80\tunmeasured\n"+
			"scalarbounds\t"+bounds+"\nscalarboundssha80\t"+boundsIdentity+"\n"+
			"top\t50\nk\t5\nseed\t0\nwidth\t1024\nddlreject\t0\nddlaccept\t1025\n"+
			"cosine\tranking\nrankingfeatures\tcombined25\nscalarscale\t0.1\nscalargates\tf\nscalarz\t0\n"+
			"identityacceptfloor\t30\nblosumacceptfloor\t0\npositional\tt\npartitions\t4\noverlap\t2\n"+
			"indexearlyreject\tf\nmarginscale\t1\n";
		final String fixedSettings=settings+profileBinding(fixedProfiles);
		final String quadrantSettings=settings+profileBinding(quadrantProfiles)+"positionmode\tfixed-quadrant\n";
		final FamilyTailCascade.LoadedClassifier omitted=load100(root,"omitted",fixedSettings);
		write(root.resolve("omitted.decision.tsv"),
			new String(decision(omitted.newWorker().assign("q",query)),StandardCharsets.US_ASCII));
		if(args.length==4){
			check(mode(omitted)==ProteinPositionMode.FIXED,"Retained parent must keep its original fixed default");
			write(root.resolve("PASS"),"STOP_TAIL_POSITION_MODE_LEGACY_PASS\n");
			System.out.println("STOP_TAIL_POSITION_MODE_LEGACY_PASS"); return;
		}
		final FamilyTailCascade.LoadedClassifier explicit=load100(root,"fixed",fixedSettings+"positionmode\tfixed\n");
		check(mode(omitted)==ProteinPositionMode.FIXED && mode(explicit)==ProteinPositionMode.FIXED,
			"Omitted selector retains fixed rather than changing the historical default");
		check(Arrays.equals(decision(omitted.newWorker().assign("q",query)),decision(explicit.newWorker().assign("q",query))),
			"Explicit fixed and omitted selector retain identical classification and work fields");
		final FamilyTailCascade.LoadedClassifier quad=load100(root,"quadrant",quadrantSettings);
		check(mode(quad)==ProteinPositionMode.FIXED_QUADRANT && quad.matchesProfiles(new char[][]{quadrant}) &&
			!quad.matchesProfiles(new char[][]{fixed}),"The caller must load actual quadrant registers, not fixed ones");
		final HashMap<String,String> directOptions=new HashMap<String,String>();
		for(String line:quadrantSettings.split("\n")){
			final String[] fields=line.split("\t",2);
			if(!fields[0].equals("marginscale")){directOptions.put(fields[0],fields[1]);}
		}
		final FamilyTailCascade.LoadedClassifier direct=FamilyTailCascade.loadCallableStructural(directOptions,1,
			new ddl.DDLFamilySearch.StructuralInputs());
		final FamilyTailCascade.Worker worker=quad.newWorker();
		check(workerMode(worker)==ProteinPositionMode.FIXED_QUADRANT,
			"The actual caller query worker must retain the selector");
		check(Arrays.equals(decision(worker.assign("q",query)),decision(direct.newWorker().assign("q",query))),
			"Caller quadrant classification must equal the direct native quadrant kernel");
		reject(()->load100(root,"unknown",fixedSettings+"positionmode\tbogus\n"),"Unknown protein position mode");
		reject(()->load100(root,"duplicate",fixedSettings+"positionmode\tfixed\npositionmode\tfixed\n"),
			"duplicate classifier binding");
		reject(()->load100(root,"empty",fixedSettings+"positionmode\t\n"),"classifier binding");
		reject(()->load100(root,"contradiction",quadrantSettings.replace("positional\tt\n","positional\tf\n")),
			"Contradictory protein position mode");
		reject(()->load100(root,"wrong-fixed-profile",fixedSettings+"positionmode\tfixed-quadrant\n"),"sketch");
		reject(()->load100(root,"wrong-quadrant-profile",
			settings+profileBinding(quadrantProfiles)+"positionmode\tfixed\n"),"sketch");
		final Path medians=write(root.resolve("medians.tsv"),
			"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\nq\t1\t100\t100\t100\t1\n");
		final Path hydro=write(root.resolve("hydro.tsv"),"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tq\t1\t0\n");
		final String ranking="heuristicarm\thitfourth_h020_c150\nfamilymedians\t"+medians+
			"\nfamilymedianssha80\tunmeasured\nfamilyhydro\t"+hydro+"\nfamilyhydrosha80\tunmeasured\n";
		for(int order=0;order<2;order++){
			final String longSettings=order==0 ? fixedSettings+"positionmode\tfixed\n" : quadrantSettings;
			final String shortSettings=order==0 ? quadrantSettings : fixedSettings+"positionmode\tfixed\n";
			final Path child=write(root.resolve("child-"+order+".tsv"),"#bbgap_tail_cascade\t2\n"+shortSettings+ranking+
				"minqueryaa\t40\nmaxqueryaa\t99\n");
			final Path parent=write(root.resolve("parent-"+order+".tsv"),"#bbgap_tail_cascade\t7\n"+longSettings+ranking+
				"minqueryaa\t100\nmaxqueryaa\t100\nshort40manifest\t"+child+"\nshort40sha80\tunmeasured\n");
			final StopTailCascadeResources composite=StopTailCascadeResources.readStructural(parent.toString(),"unmeasured");
			final FamilyTailCascade.LoadedClassifier longView=composite.load(1),shortView=composite.loadShort40(1);
			check(mode(longView)==(order==0 ? ProteinPositionMode.FIXED : ProteinPositionMode.FIXED_QUADRANT) &&
				mode(shortView)==(order==0 ? ProteinPositionMode.FIXED_QUADRANT : ProteinPositionMode.FIXED),
				"Root and SHORT selectors must remain independent through real model loading");
			check(workerMode(longView.newWorker())==mode(longView) && workerMode(shortView.newWorker())==mode(shortView),
				"Both composite query workers must retain their own selector");
		}
		final String pass="STOP_TAIL_POSITION_MODE_PASS checks="+checks
			+" max_aa=100 legacy_fixed_parity=true quadrant_profile_and_query=true independent_views=true\n";
		write(root.resolve("PASS"),pass); System.out.print(pass);
	}

	private static FamilyTailCascade.LoadedClassifier load100(Path root,String name,String settings)throws Exception{
		final Path manifest=write(root.resolve(name+"-manifest.tsv"),
			"#bbgap_tail_cascade\t2\n"+settings+"minqueryaa\t100\nmaxqueryaa\t100\n");
		return StopTailCascadeResources.readFixed100Structural(manifest.toString(),"unmeasured").load(1);
	}
	private static String profileBinding(Path list){return "ddlprofiles\t"+list+"\nddlprofilessha80\tunmeasured\n";}
	private static char[] sketch(byte[] query,ProteinKmerGeometry geometry,ProteinPositionMode mode){
		final char[] out=new char[1024]; DDLFamilyProfiles.sketchPositional(query,geometry,0,4,2,mode,out,null); return out;
	}
	private static Path profiles(Path root,String name,char[] registers,ProteinPositionMode mode)throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name+"-profiles"));
		final ByteBuilder data=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl()
			.append("tail\tq\t1024\t1\t100\t96\tnonempty\t");
		for(int kind=0;kind<3;kind++){
			if(kind>0){data.tab();}
			for(int i=0;i<registers.length;i++){
				if(i>0){data.append(',');} data.append(kind==0 ? (int)registers[i] : registers[i]==0 ? 0 : 1);
			}
		}
		final Path profile=write(dir.resolve("profile.tsv"),data.nl().toString());
		write(dir.resolve("summary.tsv"),
			"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\n"
			+"tail\t5\t0\t6\t16\t"+mode.sketch+"\t4\t2\t"+DDLFamilyProfiles.positionalTagBits(4)
			+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+"\n");
		return write(dir.resolve("list.tsv"),profile+"\tunmeasured\tunmeasured\n");
	}
	/** Inspect existing private policy state rather than add a production API just for a fixture. */
	private static ProteinPositionMode mode(FamilyTailCascade.LoadedClassifier c)throws Exception{return configMode(c);}
	private static ProteinPositionMode workerMode(FamilyTailCascade.Worker w)throws Exception{return configMode(w);}
	private static ProteinPositionMode configMode(Object object)throws Exception{
		final Field binding=object.getClass().getDeclaredField("config"); binding.setAccessible(true);
		final Object config=binding.get(object);
		final Field mode=config.getClass().getDeclaredField("positionMode");
		mode.setAccessible(true);
		return (ProteinPositionMode)mode.get(config);
	}
	private static byte[] decision(FamilyTailCascade.Assignment a){
		return (a.assigned()+"\t"+a.familyId+"\t"+a.terminal+"\t"+a.rosterRank+"\t"+a.ddlMatches+"\t"+a.ddlPosition+
			"\t"+a.dpAttempt+"\t"+a.dpAttempts+"\t"+a.score+"\t"+a.identity+"\t"+a.blosum+"\t"+
			a.scoreMeasured+"\t"+a.identityMeasured+"\t"+a.blosumMeasured+"\n").getBytes(StandardCharsets.US_ASCII);
	}
	private static Path write(Path p,String text)throws Exception{
		return Files.write(p,text.getBytes(StandardCharsets.UTF_8),java.nio.file.StandardOpenOption.CREATE_NEW);
	}
	private interface Action{void run()throws Exception;}
	private static void reject(Action action,String message)throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null
				&& expected.getMessage().toLowerCase(java.util.Locale.ROOT).contains(message.toLowerCase(java.util.Locale.ROOT)),
				"Expected rejection containing "+message+", got "+expected); return;
		}
		throw new AssertionError("Expected rejection: "+message);
	}
	private static void check(boolean okay,String why){if(!okay){throw new AssertionError(why);} checks++;}
	private static int checks;
}
