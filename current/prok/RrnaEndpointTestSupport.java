package prok;

import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;

/** Small text-fixture reader; no dependency on research comparison tools.
 * @author Raiden
 */
final class RrnaEndpointTestSupport {
	static String read(String path){
		final ByteFile in=RrnaResourceIO.open(path);final StringBuilder b=new StringBuilder();
		try{for(byte[] line;(line=in.nextLine())!=null;){b.append(new String(line,StandardCharsets.UTF_8)).append('\n');}}
		finally{if(in.close()){throw new AssertionError("Fixture read failed");}}return b.toString();
	}
	private RrnaEndpointTestSupport(){}
}
