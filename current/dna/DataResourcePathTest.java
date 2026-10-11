package dna;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPOutputStream;

import fileIO.ReadWrite;
import shared.Resources;

/** Isolated resource lookup regression; root must be a fresh disposable runtime.
 * @author Keqing
 */
public final class DataResourcePathTest {
	public static void main(String[] args) throws Exception{
		String root=null;
		boolean baseline=false, jar=false;
		for(String arg:args){
			if(arg.startsWith("root=")){root=arg.substring(5);}
			else if(arg.equals("baseline=t")){baseline=true;}
			else if(arg.equals("jar=t")){jar=true;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(root==null){throw new IllegalArgumentException("Require root=fresh_disposable_runtime");}
		final Path base=Paths.get(root).toAbsolutePath().normalize();
		require(Paths.get(Data.ROOT()).toAbsolutePath().normalize().equals(base), "Run classes from root/current so fixture writes stay isolated");
		final DataResourcePathTest test=new DataResourcePathTest(base, baseline);
		test.run(jar);
		System.out.println("RESOURCE_PATH_PASS checks="+test.checks+" baseline="+baseline+" jar="+jar);
	}

	private DataResourcePathTest(Path base_, boolean baseline_){base=base_;baseline=baseline_;}

	private void run(boolean jar) throws Exception{
		final Path resources=base.resolve("resources"), networks=base.resolve("networks");
		final Path dori=base.resolve("dori/resources"), doriNets=base.resolve("dori/networks");
		final Path nersc=base.resolve("nersc/resources"), nerscNets=base.resolve("nersc/networks");
		Data.DORI_RESOURCES=dori+"/"; Data.DORI_NETWORKS=doriNets+"/";
		Data.NERSC_RESOURCES=nersc+"/"; Data.NERSC_NETWORKS=nerscNets+"/";
		final PrintStream prior=System.err;
		final ByteArrayOutputStream warnings=new ByteArrayOutputStream();
		try{
			System.setErr(new PrintStream(warnings, true, "UTF-8"));
			final Path exact=put(resources,"audit-exact.txt",false);
			put(resources,"audit-exact.txt.gz",true);
			equal(Data.findPath("?audit-exact.txt",true),exact,"Exact name wins within a location");
			final Path gz=put(resources,"audit-gzip.txt.gz",true);
			equal(Data.findPath("?audit-gzip.txt",true),baseline?null:gz,"Compressed local resource");
			equal(Data.findPath("?audit-gzip.txt.gz",false),gz,"Explicit gzip name");
			put(resources,"audit-double.txt.gz.gz",true);
			equal(Data.findPath("?audit-double.txt.gz",false),null,"No double gzip suffix");
			put(resources,"audit-double.txt.gzip.gz",true);
			equal(Data.findPath("?audit-double.txt.gzip",false),null,"No suffix after gzip extension");
			put(resources,"audit-upper.GZ.gz",true);
			equal(Data.findPath("?audit-upper.GZ",false),null,"Compressed-extension check ignores case");
			equal(Data.findPath("literal-missing.txt",true),Paths.get("literal-missing.txt"),"Literal path is returned unchanged");
			equal(Data.findPath(gz.toString().substring(0,gz.toString().length()-3),true),resources.resolve("audit-gzip.txt"),"Literal path does not gain suffix");
			final Path nested=put(resources,"nested/audit space.txt.gz",true);
			equal(Data.findPath("?nested/audit space.txt",false),baseline?null:nested,"Nested resource with space");
			final Path local=put(resources,"audit-priority.txt.gz",true);
			final Path shared=put(dori,"audit-priority.txt",false);
			equal(Data.findPath("?audit-priority.txt",false),baseline?shared:local,"Local compressed beats shared exact");
			final Path rootGz=put(base,"audit-root.txt.gz",true);
			equal(Data.findPath("?audit-root.txt",false),baseline?null:rootGz,"Root fallback accepts gzip");
			final Path network=put(networks,"audit-net.bbnet.gz",true);
			final Path compatibility=put(resources,"audit-net.bbnet",false);
			equal(Data.findPath("?audit-net.bbnet",false),baseline?compatibility:network,"Networks precede compatibility resources");
			final Path networkExact=put(networks,"audit-net.bbnet",false);
			equal(Data.findPath("?audit-net.bbnet",false),networkExact,"Exact network wins over compressed network");
			final Path bundle=put(resources,"audit-bundle.bbnets.gz",true);
			equal(Data.findPath("?audit-bundle.bbnets",false),baseline?null:bundle,"Network bundle compatibility directory");
			final Path dr=put(dori,"audit-dori.txt.gz",true);
			equal(Data.findPath("?audit-dori.txt",false),baseline?null:dr,"Shared Dori resource gzip");
			final Path dn=put(doriNets,"audit-dori.bbnet.gz",true);
			equal(Data.findPath("?audit-dori.bbnet",false),baseline?null:dn,"Shared Dori network gzip");
			final Path dc=put(dori,"audit-compat.bbnet.gz",true);
			equal(Data.findPath("?audit-compat.bbnet",false),baseline?null:dc,"Shared Dori compatibility gzip");
			final Path nr=put(nersc,"audit-nersc.txt.gz",true);
			equal(Data.findPath("?audit-nersc.txt",false),baseline?null:nr,"Shared NERSC resource gzip");
			final Path nn=put(nerscNets,"audit-nersc.bbnet.gz",true);
			equal(Data.findPath("?audit-nersc.bbnet",false),baseline?null:nn,"Shared NERSC network gzip");
			final Path nc=put(nersc,"audit-nersc-compat.bbnet.gz",true);
			equal(Data.findPath("?audit-nersc-compat.bbnet",false),baseline?null:nc,"Shared NERSC compatibility gzip");
			if(!baseline){
				equal(Resources.find("?audit-gzip.txt"),gz,"Shared Resources wrapper inherits fallback");
				require(ReadWrite.readString(Data.findPath("?audit-gzip.txt",false)).equals(CONTENT),"Resolved gzip must decode through BBTools I/O");checks++;
			}
			if(jar){
				final String classpath=Data.findPath("?audit-classpath.txt",false);
				require(baseline?classpath==null:classpath!=null && classpath.startsWith("jar:") && classpath.endsWith("!/audit-classpath.txt.gz"),"JAR resource gzip fallback: "+classpath);checks++;
				final String netpath=Data.findPath("?audit-classpath.bbnet",false);
				require(baseline?netpath==null:netpath!=null && netpath.startsWith("jar:") && netpath.endsWith("!/networks/audit-classpath.bbnet.gz"),"JAR network gzip fallback: "+netpath);checks++;
			}
			if(!baseline){require(warnings.size()==0,"Successful fallback must not emit a missing-resource warning: "+warnings.toString("UTF-8"));checks++;}
			warnings.reset();
			equal(Data.findPath("?audit-absent.txt",false),null,"Quiet missing resource");
			require(warnings.size()==0,"warn=false must remain quiet");checks++;
			equal(Data.findPath("?audit-absent.txt",true),null,"Warning missing resource");
			final String warning=warnings.toString("UTF-8"), marker="Warning!  Cannot find";
			require(warning.contains(marker) && warning.indexOf(marker)==warning.lastIndexOf(marker),"Exactly one missing-resource warning expected");checks++;
		}finally{System.setErr(prior);}
	}

	private Path put(Path directory,String name,boolean gzip) throws IOException{
		final Path target=directory.resolve(name);
		if(Files.exists(target)){throw new IOException("Fixture already exists: "+target);}
		Files.createDirectories(target.getParent());
		final byte[] bytes=CONTENT.getBytes(StandardCharsets.UTF_8);
		if(gzip){try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(target))){out.write(bytes);}}
		else{Files.write(target,bytes);}
		return target;
	}
	private void equal(String actual,Path expected,String why){
		require(expected==null?actual==null:expected.toString().equals(actual),why+": expected="+expected+", actual="+actual);checks++;
	}
	private static void require(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	private static final String CONTENT="resource gzip fixture\n";
	private final Path base;
	private final boolean baseline;
	private int checks;
}
