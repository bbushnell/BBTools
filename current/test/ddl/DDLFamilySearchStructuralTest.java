package test.ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLFamilySearch.Index;
import ddl.DDLFamilySearch.StructuralInputs;
import ddl.DDLFamilySearch.Workspace;
import ddl.ProteinKmerGeometry;

/** Structural loading accepts inherited pins or explicit unmeasured markers while guarding files.
 * @author Yelan
 */
public final class DDLFamilySearchStructuralTest{

	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("ddl_structural_");
		final Path good=fixture(root.resolve("good"),summary("",""),profiles(4));
		final Path firstList=list(root.resolve("first.tsv"),good,PIN_A,PIN_B);
		final Path secondList=list(root.resolve("second.tsv"),good,PIN_B,PIN_A);

		final StructuralInputs firstInputs=new StructuralInputs(), secondInputs=new StructuralInputs();
		final Index expected=DDLFamilySearch.loadForAttributionStructural(firstList,4,ProteinKmerGeometry.exact(5),0,
			DDLFamilyProfiles.SKETCH_PLAIN,4,8,firstInputs);
		final Index structural=DDLFamilySearch.loadForAttributionStructural(secondList,4,ProteinKmerGeometry.exact(5),0,
			DDLFamilyProfiles.SKETCH_PLAIN,4,8,secondInputs);
		firstInputs.verify(); secondInputs.verify();
		check(Arrays.equals(expected.names,structural.names),"Structural loader preserves roster order");
		for(int i=0; i<expected.families(); i++){
			check(Arrays.equals(expected.profileCopy(i),structural.profileCopy(i)),"Structural rows equal inherited-pin rows");
		}
		final Workspace a=new Workspace(expected.families()), b=new Workspace(structural.families());
		final char[] query={(char)1,(char)2,(char)3,(char)4};
		expected.score(query,a); structural.score(query,b);
		final int[] rankA=new int[4], rankB=new int[4];
		check(expected.top(a,rankA)==structural.top(b,rankB) && Arrays.equals(rankA,rankB),
			"Structural scoring and ranking match inherited-pin loading");

		final String[][] identities={{"unmeasured","unmeasured"},{"unmeasured",PIN_B},{PIN_A,"unmeasured"}};
		for(int i=0;i<identities.length;i++){
			final StructuralInputs inputs=new StructuralInputs();
			final Index loaded=DDLFamilySearch.loadForAttributionStructural(
				list(root.resolve("unmeasured"+i+".tsv"),good,identities[i][0],identities[i][1]),
				4,ProteinKmerGeometry.exact(5),0,DDLFamilyProfiles.SKETCH_PLAIN,4,8,inputs);
			inputs.verify();
			check(Arrays.equals(expected.names,loaded.names),"Unmeasured identity preserves roster order");
			for(int family=0;family<expected.families();family++){
				check(Arrays.equals(expected.profileCopy(family),loaded.profileCopy(family)),"Unmeasured identity preserves rows");
			}
			final Workspace work=new Workspace(loaded.families());
			loaded.score(query,work);
			check(Arrays.equals(a.scores,work.scores),"Unmeasured identity preserves scores");
			check(expected.top(a,rankA)==loaded.top(work,rankB) && Arrays.equals(rankA,rankB),
				"Unmeasured identity preserves ranking");
		}
		int invalidCase=0;
		for(final String invalid:new String[]{"","UNMEASURED","unmeasured "," unmeasured","unmeasuredx","null","1111111111111111111"}){
			for(int column=0;column<2;column++){
				final Path invalidList=list(root.resolve("invalid"+(invalidCase++)+".tsv"),good,
					column==0 ? invalid : "unmeasured",column==1 ? invalid : "unmeasured");
				expectFailure(new Action(){@Override public void run() throws Exception{
					DDLFamilySearch.loadForAttributionStructural(invalidList,4,ProteinKmerGeometry.exact(5),0,
						DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
				}},"Only literal unmeasured or a complete inherited pin is allowed in either identity field");
			}
		}

		expectFailure(new Action(){@Override public void run() throws Exception{
			DDLFamilySearch.loadForAttributionStructural(list(root.resolve("badpin.tsv"),good,"bad",PIN_B),
				4,ProteinKmerGeometry.exact(5),0,DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Structural loading still validates declared pin syntax");
		expectFailure(new Action(){@Override public void run() throws Exception{
			final Path malformed=fixture(root.resolve("malformed"),"bad\n",profiles(4));
			DDLFamilySearch.loadForAttributionStructural(list(root.resolve("malformed.tsv"),malformed,
				"unmeasured","unmeasured"),4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Malformed profile summaries are still rejected");
		expectFailure(new Action(){@Override public void run() throws Exception{
			final Path mismatch=fixture(root.resolve("geometry"),summary("\tgeometry","\texact6"),profiles(4));
			DDLFamilySearch.loadForAttributionStructural(list(root.resolve("geometry.tsv"),mismatch,
				"unmeasured","unmeasured"),4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Geometry metadata must still match the caller");
		expectFailure(new Action(){@Override public void run() throws Exception{
			Files.write(root.resolve("duplicate.tsv"),(good+"\t"+PIN_A+"\t"+PIN_B+"\n"+
				good+"\t"+PIN_A+"\t"+PIN_B+"\n").getBytes(StandardCharsets.US_ASCII));
			DDLFamilySearch.loadForAttributionStructural(root.resolve("duplicate.tsv"),4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Duplicate profile-list rows remain invalid");
		expectFailure(new Action(){@Override public void run() throws Exception{
			DDLFamilySearch.loadForAttributionStructural(list(root.resolve("missing.tsv"),root.resolve("missing/profiles.tsv"),
				PIN_A,PIN_B),4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Missing nested input files are rejected before reads");
		expectFailure(new Action(){@Override public void run() throws Exception{
			final Path duplicateFamily=fixture(root.resolve("duplicate_family"),summary("",""),
				profiles(4)+row("a",new int[]{1,2,0,4}));
			DDLFamilySearch.loadForAttributionStructural(list(root.resolve("duplicate_family.tsv"),duplicateFamily,
				"unmeasured","unmeasured"),4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,new StructuralInputs());
		}},"Unmeasured identities do not permit duplicate family rows");
		for(int component=0;component<3;component++){
			final Path profile=fixture(root.resolve("change"+component),summary("",""),profiles(4));
			final Path profileList=list(root.resolve("change"+component+".tsv"),profile,"unmeasured","unmeasured");
			final StructuralInputs inputs=new StructuralInputs();
			DDLFamilySearch.loadForAttributionStructural(profileList,4,ProteinKmerGeometry.exact(5),0,
				DDLFamilyProfiles.SKETCH_PLAIN,4,8,inputs);
			final Path changedFile=component==0 ? profileList : component==1 ? profile : profile.getParent().resolve("summary.tsv");
			Files.write(changedFile,"#changed\n".getBytes(StandardCharsets.US_ASCII),java.nio.file.StandardOpenOption.APPEND);
			expectFailure(new Action(){@Override public void run() throws Exception{inputs.verify();}},
				"Unmeasured list, profile and summary remain tracked after loading");
		}

		final StructuralInputs changed=new StructuralInputs();
		changed.track(good);
		Files.write(good,(profiles(4)+"#changed\n").getBytes(StandardCharsets.US_ASCII));
		expectFailure(new Action(){@Override public void run() throws Exception{changed.verify();}},
			"Tracked metadata changes must be detected");
		expectFailure(new Action(){@Override public void run() throws Exception{changed.track(good);}},
			"Repeated track calls verify rather than replacing the baseline");

		System.out.println("DDL_FAMILY_SEARCH_STRUCTURAL_TEST_PASS");
	}

	private static Path fixture(Path dir,String summary,String profiles) throws Exception{
		Files.createDirectories(dir);
		Files.write(dir.resolve("summary.tsv"),summary.getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("profiles.tsv"),profiles.getBytes(StandardCharsets.US_ASCII));
		return dir.resolve("profiles.tsv");
	}

	private static Path list(Path list,Path profiles,String profilePin,String summaryPin) throws Exception{
		final Path parent=list.getParent();
		if(parent!=null){Files.createDirectories(parent);}
		Files.write(list,(profiles+"\t"+profilePin+"\t"+summaryPin+"\n").getBytes(StandardCharsets.US_ASCII));
		return list;
	}

	private static String summary(String extraHeader,String extraValue){
		return "namespace\tk\tseed\texponent"+extraHeader+"\nfixture\t5\t0\t6"+extraValue+"\n";
	}

	private static String profiles(int width){
		return DDLFamilyProfiles.HEADER16+"\n"+row("a",new int[]{1,2,0,4})+row("b",new int[]{1,0,3,0});
	}

	private static String row(String id,int[] values){
		final StringBuilder line=new StringBuilder("fixture\t").append(id).append('\t').append(values.length)
			.append("\t2\t16\t12\tnonempty");
		for(int column=0; column<3; column++){
			line.append('\t');
			for(int i=0; i<values.length; i++){
				if(i>0){line.append(',');}
				line.append(column==0 ? values[i] : values[i]==0 ? 0 : 1);
			}
		}
		return line.append('\n').toString();
	}

	private static void expectFailure(Action action,String message) throws Exception{
		boolean failed=false;
		try{action.run();}
		catch(RuntimeException expected){failed=true;}
		catch(java.io.IOException expected){failed=true;}
		check(failed,message);
	}

	private static void check(boolean valid,String message){if(!valid){throw new AssertionError(message);}}

	private interface Action{
		void run() throws Exception;
	}

	private static final String PIN_A="11111111111111111111";
	private static final String PIN_B="22222222222222222222";
}
