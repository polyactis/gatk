/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk;

import java.util.*;
import java.io.InputStream;
import java.io.IOException;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.filters.FilterManager;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.PluginManager;
import org.broadinstitute.sting.utils.help.DisplayNameTaglet;
import org.apache.log4j.Logger;
import net.sf.picard.filter.SamRecordFilter;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 17, 2009
 * Time: 3:14:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class WalkerManager extends PluginManager<Walker> {

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(WalkerManager.class);

    /**
     * A collection of help text for walkers and their enclosing packages.
     */
    private Properties helpText = new Properties();

    public WalkerManager() {
        super(Walker.class,"walker","Walker");
        InputStream helpSourceFile = getClass().getClassLoader().getResourceAsStream("help.properties");
        if(helpSourceFile != null) {
            try {
                helpText.load(helpSourceFile);
            }
            catch(IOException ex) {
                throw new StingException("Unable to process help data");
            }
        }
    }

    /**
     * Get the list of walkers currently available to the GATK, organized
     * by package.
     * @return Names of currently available walkers.
     */
    public Map<String,Collection<Class<? extends Walker>>> getWalkerNamesByPackage() {
        Map<String,Collection<Class<? extends Walker>>> walkersByPackage = new HashMap<String,Collection<Class<? extends Walker>>>();
        for(Class<? extends Walker> walker: pluginsByName.values()) {
            String walkerPackage = walker.getPackage() != null ? walker.getPackage().getName() : "<unpackaged>";
            if(!walkersByPackage.containsKey(walkerPackage))
                walkersByPackage.put(walkerPackage,new ArrayList<Class<? extends Walker>>());
            walkersByPackage.get(walkerPackage).add(walker);
        }
        return Collections.unmodifiableMap(walkersByPackage);
    }

    /**
     * Gets the display name for a given package.
     * @param packageName Fully qualified package name.
     * @return A suitable display name for the package.
     */
    public String getPackageDisplayName(String packageName) {
        String specifiedDisplayName = helpText.getProperty(packageName+"."+ DisplayNameTaglet.NAME);
        return specifiedDisplayName != null ? specifiedDisplayName : packageName.substring(packageName.lastIndexOf('.')+1);
    }

    /**
     * Gets the help text associated with a given package name.
     * @param packageName Package for which to search for help text.
     * @return Package help text, or "" if none exists.
     */
    public String getPackageHelpText(String packageName) {
        if(!helpText.containsKey(packageName))
            return "";
        return helpText.getProperty(packageName);    
    }

    /**
     * Gets the help text associated with a given walker type.
     * @param walkerType Type of walker for which to search for help text.
     * @return Package help text, or "" if none exists.
     */
    public String getWalkerHelpText(Class<? extends Walker> walkerType) {
        String walkerName = walkerType.getName();
        if(!helpText.containsKey(walkerName))
            return "";
        return helpText.getProperty(walkerName);
    }

    /**
     * Retrieves the walker class given a walker name.
     * @param walkerName Name of the walker.
     * @return Class representing the walker.
     */
    public Class getWalkerClassByName(String walkerName) {
        return pluginsByName.get(walkerName);
    }

    /**
     * Gets the data source for the provided walker.
     * @param walker The walker.
     * @return Which type of data source to traverse over...reads or reference?
     */
    public static DataSource getWalkerDataSource(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        By byDataSource = walkerClass.getAnnotation(By.class);
        if( byDataSource == null )
            throw new StingException("Unable to find By annotation for walker class " + walkerClass.getName());
        return byDataSource.value();
    }

    /**
     * Determine whether the given walker supports the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for .
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, DataSource dataSource) {
        Allows allowsDataSource = getWalkerAllowed(walker);

        // Allows is less restrictive than requires.  If an allows
        // clause is not specified, any kind of data is allowed.
        if( allowsDataSource == null )
            return true;

        return Arrays.asList(allowsDataSource.value()).contains(dataSource);
    }

    /**
     * Determine whether the given walker supports the given reference ordered data.
     * @param walker Walker to query.
     * @param rod Source to check.
     * @return True if the walker forbids this data type.  False otherwise.
     */
    public static boolean isAllowed(Walker walker, ReferenceOrderedData<? extends ReferenceOrderedDatum> rod) {
        Allows allowsDataSource = getWalkerAllowed(walker);

        // Allows is less restrictive than requires.  If an allows
        // clause is not specified, any kind of data is allowed.
        if( allowsDataSource == null )
            return true;

        // The difference between unspecified RMD and the empty set of metadata can't be detected.
        // Treat an empty 'allows' as 'allow everything'.  Maybe we can have a special RMD flag to account for this
        // case in the future.
        if( allowsDataSource.referenceMetaData().length == 0 )
            return true;

        for( RMD allowed: allowsDataSource.referenceMetaData() ) {
            if( rod.matches(allowed.name(),allowed.type()) )
                return true;
        }
        return false;
    }

    /**
     * Determine whether the given walker requires the given data source.
     * @param walker Walker to query.
     * @param dataSource Source to check for.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static boolean isRequired(Walker walker, DataSource dataSource) {
        Requires requiresDataSource = getWalkerRequirements(walker);
        return Arrays.asList(requiresDataSource.value()).contains(dataSource);
    }

    /**
     * Get a list of RODs required by the walker.
     * @param walker Walker to query.
     * @return True if the walker allows this data type.  False otherwise.
     */
    public static List<RMD> getRequiredMetaData(Walker walker) {
        Requires requiresDataSource = getWalkerRequirements(walker);
        return Arrays.asList(requiresDataSource.referenceMetaData());
    }

    /**
     * Extracts filters that the walker has requested be run on the dataset.
     * @param walker Walker to inspect for filtering requests.
     * @param filterManager Manages the creation of filters.
     * @return A non-empty list of filters to apply to the reads.
     */
    public static List<SamRecordFilter> getReadFilters(Walker walker, FilterManager filterManager) {
        List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
        for(Class<? extends SamRecordFilter> filterType: getReadFilterTypes(walker))
            filters.add(filterManager.createFilterByType(filterType));
        return filters;
    }

    /**
     * Create a name for this type of walker.
     *
     * @param walkerType The type of walker.
     * @return A name for this type of walker.
     */
    @Override
    public String getName(Class<? extends Walker> walkerType) {
        String walkerName = "";

        if (walkerType.getAnnotation(WalkerName.class) != null)
            walkerName = walkerType.getAnnotation(WalkerName.class).value().trim();
        else
            walkerName = super.getName(walkerType);

        return walkerName;
    }

    /**
     * Utility to get the requires attribute from the walker.
     * Throws an exception if requirements are missing.
     * @param walker Walker to query for required data.
     * @return Required data attribute.
     */
    private static Requires getWalkerRequirements(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        Requires requiresDataSource = walkerClass.getAnnotation(Requires.class);
        if( requiresDataSource == null )
            throw new StingException( "Unable to find data types required by walker class " + walkerClass.getName());
        return requiresDataSource;
    }

    /**
     * Utility to get the forbidden attribute from the walker.
     * @param walker Walker to query for required data.
     * @return Required data attribute.  Null if forbidden info isn't present.
     */
    private static Allows getWalkerAllowed(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        Allows allowsDataSource = walkerClass.getAnnotation(Allows.class);
        return allowsDataSource;
    }

    /**
     * Gets the list of filtering classes specified as walker annotations.
     * @param walker The walker to inspect.
     * @return An array of types extending from SamRecordFilter.  Will never be null.
     */
    private static Class<? extends SamRecordFilter>[] getReadFilterTypes(Walker walker) {
        Class<? extends Walker> walkerClass = walker.getClass();
        if( !walkerClass.isAnnotationPresent(ReadFilters.class) )
            return new Class[0];
        return walkerClass.getAnnotation(ReadFilters.class).value();
    }
}
