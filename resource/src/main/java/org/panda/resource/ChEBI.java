package org.panda.resource;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Provides ChEBI ID to name conversion.
 *
 * @author Ozgun Babur
 */
public class ChEBI extends FileServer
{
	private Map<String, String> idToName;

	private static ChEBI instance;

	public static ChEBI get()
	{
		if (instance == null) instance = new ChEBI();
		return instance;
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{
			"ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz",
			"ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"};
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"ChEBI-names1.txt", "ChEBI-names2.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		idToName = new HashMap<>();
		getResourceAsStream(getLocalFilenames()[0]).skip(1).forEach(line ->
		{
			String[] token = line.split("\t");

			if (token[2].equals("NAME"))
			{
				idToName.put("CHEBI:" + token[1], token[4]);
			}
		});

		getResourceAsStream(getLocalFilenames()[1]).skip(1).forEach(line ->
		{
			String[] token = line.split("\t");

			if (token.length > 5 && !token[5].isEmpty() && !token[5].equals("null"))
			{
				idToName.put(token[2], token[5]);
			}
		});
		return true;
	}

	/**
	 * @param id ChEBI ID that starts with "CHEBI:"
	 */
	public String getName(String id)
	{
		return idToName.get(id);
	}

	public Map<String, String> getIdToNameMapping()
	{
		return idToName;
	}

	public static void main(String[] args) throws IOException
	{
//		for (String s : "CHEBI:17295, CHEBI:16977, CHEBI:17196, CHEBI:17191, CHEBI:16414, CHEBI:17203, CHEBI:16828, CHEBI:15603, CHEBI:16015".split(", "))
//		{
//			System.out.println(s + "\t" + get().getName(s));
//		}

		replaceIDsWithNameInSIF("/home/babur/Documents/Analyses/MrOS/links-between-metabolites-and-proteins.sif");
	}

	private static void replaceIDsWithNameInSIF(String filename) throws IOException
	{
		List<String> lines = new ArrayList<>();
		Files.lines(Paths.get(filename)).map(l -> l.split("\t")).forEach(t ->
		{
			if (t[0].startsWith("CHEBI:")) t[0] = get().getName(t[0]) != null ? get().getName(t[0]) : t[0];
			if (t[2].startsWith("CHEBI:")) t[2] = get().getName(t[2]) != null ? get().getName(t[2]) : t[2];
			String s = t[0];
			for (int i = 1; i < t.length; i++)
			{
				s += "\t" + t[i];
			}
			lines.add(s);
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		for (String line : lines)
		{
			writer.write(line + "\n");
		}
		writer.close();
	}
}
