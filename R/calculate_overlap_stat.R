


#' Get Matched Matrix
#'
#' Converts a full matched table into a matrix where each row
#' corresponds to a treated unit and each column contains the IDs of
#' matched control units. Missing values (if any) are filled with
#' `NA`. Fixed version that properly handles control reuse across
#' subclasses.
#'
#' @param full_matched_table A data frame containing the matched data.
#'   Must include the columns:
#' \itemize{
#'   \item{\code{id}: Unique identifiers for each unit.}
#'   \item{\code{Z}: Indicator for treated (\code{1}) or control (\code{0}) units.}
#'   \item{\code{subclass}: Subclass assignments for matching.}
#' }
#'
#' @return A matrix where: \item{Rows}{Represent treated units
#'   (indexed by their IDs).} \item{Columns}{Contain IDs of matched
#'   control units for each treated unit. Missing matches are filled
#'   with \code{NA}.}
#'
#' @export
get_matched_matrix <- function(full_matched_table) {
  # Get unique treated units
  treated_units <- unique(full_matched_table$id[full_matched_table$Z == 1])

  matched_control_ids <- list()

  # For each treated unit, get ALL controls it's matched to (across all subclasses)
  for (treated_id in treated_units) {
    # Find all subclasses this treated unit appears in
    treated_subclasses <- unique(full_matched_table$subclass[full_matched_table$id == treated_id & full_matched_table$Z == 1])

    # Get all controls from those subclasses
    all_controls <- c()
    for (subclass in treated_subclasses) {
      controls_in_subclass <- full_matched_table$id[full_matched_table$subclass == subclass & full_matched_table$Z == 0]
      all_controls <- c(all_controls, controls_in_subclass)
    }

    # Remove duplicates - this is the key fix!
    # Convert to numeric if possible for proper handling
    if (all(grepl("^[0-9]+$", all_controls))) {
      unique_controls <- unique(as.numeric(all_controls))
    } else {
      unique_controls <- unique(all_controls)
    }

    matched_control_ids[[as.character(treated_id)]] <- unique_controls
  }

  # Create matrix
  max_controls <- max(sapply(matched_control_ids, length))
  matched_matrix <- do.call(rbind, lapply(matched_control_ids, function(ids) {
    c(ids, rep(NA, max_controls - length(ids)))
  }))

  # Ensure proper row names (treated unit IDs)
  rownames(matched_matrix) <- names(matched_control_ids)

  return(matched_matrix)
}

# ============================================================================
# METRIC 1: PAIRWISE SHARED CONTROLS (Original metric, renamed for clarity)
# ============================================================================

#' Compute Pairwise Shared Controls
#'
#' Computes the number of shared control neighbors between each pair
#' of treated units based on a given matched matrix. This creates a
#' pairwise comparison matrix.
#'
#' @param matched_matrix A matrix where rows are treated units and
#'   columns contain control IDs
#'
#' @return A symmetric matrix where entry i,j is the count of shared
#'   controls between treated units i and j
#'
#' @export
compute_pairwise_shared_controls <- function(matched_matrix) {
  N1 <- nrow(matched_matrix)

  # Initialize with zeros
  shared_neighbors <- matrix(0, nrow = N1, ncol = N1)
  rownames(shared_neighbors) <- rownames(matched_matrix)
  colnames(shared_neighbors) <- rownames(matched_matrix)

  for (i in 1:(N1 - 1)) {
    for (j in (i + 1):N1) {
      # Get non-NA controls for each treated unit
      controls_i <- matched_matrix[i, ]
      controls_j <- matched_matrix[j, ]

      # Remove NAs before computing intersection
      controls_i <- controls_i[!is.na(controls_i)]
      controls_j <- controls_j[!is.na(controls_j)]

      # Find intersection
      shared <- intersect(controls_i, controls_j)
      shared_count <- length(shared)

      shared_neighbors[i, j] <- shared_neighbors[j, i] <- shared_count
    }
  }

  return(shared_neighbors)
}

#' Compute Pairwise Overlap Statistics
#'
#' Computes summary statistics for pairwise overlap of shared controls.
#' This represents how many controls are shared between pairs of treated units.
#'
#' @param pairwise_shared_matrix A matrix from compute_pairwise_shared_controls
#'
#' @return A list with median, 75th percentile, and max statistics
#'
#' @export
compute_pairwise_overlap_statistics <- function(pairwise_shared_matrix) {
  N1 <- nrow(pairwise_shared_matrix)

  shared_neighbors_binary <- pairwise_shared_matrix > 0
  n_shared_treated_vec <- rowSums(shared_neighbors_binary)
  n_shared_controls_vec <- rowSums(pairwise_shared_matrix)

  list(
    avg_shared_controls = median(n_shared_controls_vec),
    p75_shared_controls = quantile(n_shared_controls_vec, probs = 0.75),
    max_shared_controls = max(n_shared_controls_vec),
    avg_shared_treated = median(n_shared_treated_vec),
    p75_shared_treated = quantile(n_shared_treated_vec, probs = 0.75),
    max_shared_treated = max(n_shared_treated_vec)
  )
}

# ============================================================================
# METRIC 2: MEAN CONTROL REUSE (Mean of K(j))
# ============================================================================

#' Compute Mean Control Reuse
#'
#' Calculates the mean of K(j), where K(j) is the number of times each
#' control is matched to a treated unit. This metric ranges from 1 (no
#' reuse) to n_T (all controls used by all treated).
#'
#' @param matched_matrix A matrix where rows are treated units and
#'   columns contain control IDs
#'
#' @return A list containing: \item{mean_reuse}{Mean number of times a
#'   control is matched} \item{median_reuse}{Median number of times a
#'   control is matched} \item{max_reuse}{Maximum number of times a
#'   control is matched} \item{control_reuse_counts}{Named vector of
#'   reuse counts for each control}
#'
#' @export
compute_mean_control_reuse <- function(matched_matrix) {
  # Flatten the matrix and remove NAs to get all control assignments
  all_control_assignments <- as.vector(matched_matrix)
  all_control_assignments <- all_control_assignments[!is.na(all_control_assignments)]

  if (length(all_control_assignments) == 0) {
    return(list(
      mean_reuse = NA,
      median_reuse = NA,
      max_reuse = NA,
      control_reuse_counts = numeric(0)
    ))
  }

  # Count how many times each control appears
  control_counts <- table(all_control_assignments)

  # Calculate statistics
  list(
    mean_reuse = mean(control_counts),
    median_reuse = median(control_counts),
    max_reuse = max(control_counts),
    control_reuse_counts = as.numeric(control_counts)
  )
}

# ============================================================================
# METRIC 3: SHARED CONTROLS PER TREATED UNIT
# ============================================================================

#' Compute Shared Controls Per Treated
#'
#' For each treated unit, calculates how many of its matched controls
#' are shared with at least one other treated unit. This metric ranges
#' from 0 (no sharing) to k (all k matched controls are shared), where
#' k is the number of matches per treated unit.
#'
#' @param matched_matrix A matrix where rows are treated units and
#'   columns contain control IDs
#'
#' @return A list containing: \item{mean_shared}{Mean number of shared
#'   controls per treated unit} \item{median_shared}{Median number of
#'   shared controls per treated unit} \item{prop_shared}{Mean
#'   proportion of controls that are shared} \item{max_shared}{Maximum
#'   number of shared controls for any treated unit}
#'   \item{shared_per_treated}{Vector of shared control counts for
#'   each treated unit}
#'
#' @export
compute_shared_controls_per_treated <- function(matched_matrix) {
  N_treated <- nrow(matched_matrix)

  if (N_treated <= 1) {
    return(list(
      mean_shared = 0,
      median_shared = 0,
      prop_shared = 0,
      max_shared = 0,
      shared_per_treated = numeric(N_treated)
    ))
  }

  shared_counts <- numeric(N_treated)
  proportions_shared <- numeric(N_treated)

  for (i in 1:N_treated) {
    # Get controls for treated unit i
    controls_i <- matched_matrix[i, ]
    controls_i <- controls_i[!is.na(controls_i)]

    if (length(controls_i) == 0) {
      shared_counts[i] <- 0
      proportions_shared[i] <- 0
      next
    }

    # Get all controls for other treated units
    other_controls <- c()
    for (j in setdiff(1:N_treated, i)) {
      controls_j <- matched_matrix[j, ]
      controls_j <- controls_j[!is.na(controls_j)]
      other_controls <- c(other_controls, controls_j)
    }

    # Find unique controls used by others
    other_controls <- unique(other_controls)

    # Calculate intersection
    shared <- intersect(controls_i, other_controls)
    shared_counts[i] <- length(shared)
    proportions_shared[i] <- length(shared) / length(controls_i)
  }

  list(
    mean_shared = mean(shared_counts),
    median_shared = median(shared_counts),
    prop_shared = mean(proportions_shared),
    max_shared = max(shared_counts),
    shared_per_treated = shared_counts
  )
}

# ============================================================================
# MAIN FUNCTIONS
# ============================================================================

#' Calculate Overlap Statistics from Matched Table
#'
#' Calculates all three overlap metrics for a matched dataset from a
#' pre-existing full matched table.
#'
#' @param full_matched_table A data frame containing the matched data.
#'
#' @return A list containing all three sets of overlap statistics.
#'
#' @export
calculate_overlap_stats_from_table <- function(full_matched_table) {
  if (is.null(full_matched_table) || nrow(full_matched_table) == 0) {
    return(list(
      pairwise_overlap = list(
        avg_shared_controls = NA,
        p75_shared_controls = NA,
        max_shared_controls = NA,
        avg_shared_treated = NA,
        p75_shared_treated = NA,
        max_shared_treated = NA
      ),
      control_reuse = list(
        mean_reuse = NA,
        median_reuse = NA,
        max_reuse = NA
      ),
      shared_per_treated = list(
        mean_shared = NA,
        median_shared = NA,
        prop_shared = NA,
        max_shared = NA
      )
    ))
  }

  # Step 1: Get matched matrix
  matched_matrix <- tryCatch({
    get_matched_matrix(full_matched_table)
  }, error = function(e) {
    warning("Failed to get matched matrix for overlap statistics: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(matched_matrix)) {
    return(list(
      pairwise_overlap = list(
        avg_shared_controls = NA,
        p75_shared_controls = NA,
        max_shared_controls = NA,
        avg_shared_treated = NA,
        p75_shared_treated = NA,
        max_shared_treated = NA
      ),
      control_reuse = list(
        mean_reuse = NA,
        median_reuse = NA,
        max_reuse = NA
      ),
      shared_per_treated = list(
        mean_shared = NA,
        median_shared = NA,
        prop_shared = NA,
        max_shared = NA
      )
    ))
  }

  # Calculate all three metrics
  results <- list()

  # Metric 1: Pairwise overlap (original metric)
  pairwise_shared <- tryCatch({
    compute_pairwise_shared_controls(matched_matrix)
  }, error = function(e) {
    warning("Failed to compute pairwise shared controls: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (!is.null(pairwise_shared)) {
    results$pairwise_overlap <- tryCatch({
      compute_pairwise_overlap_statistics(pairwise_shared)
    }, error = function(e) {
      list(
        avg_shared_controls = NA,
        p75_shared_controls = NA,
        max_shared_controls = NA,
        avg_shared_treated = NA,
        p75_shared_treated = NA,
        max_shared_treated = NA
      )
    })
  } else {
    results$pairwise_overlap <- list(
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      max_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA,
      max_shared_treated = NA
    )
  }

  # Metric 2: Mean control reuse
  results$control_reuse <- tryCatch({
    stats <- compute_mean_control_reuse(matched_matrix)
    list(
      mean_reuse = stats$mean_reuse,
      median_reuse = stats$median_reuse,
      max_reuse = stats$max_reuse
    )
  }, error = function(e) {
    warning("Failed to compute mean control reuse: ", e$message, call. = FALSE)
    list(mean_reuse = NA, median_reuse = NA, max_reuse = NA)
  })

  # Metric 3: Shared controls per treated
  results$shared_per_treated <- tryCatch({
    stats <- compute_shared_controls_per_treated(matched_matrix)
    list(
      mean_shared = stats$mean_shared,
      median_shared = stats$median_shared,
      prop_shared = stats$prop_shared,
      max_shared = stats$max_shared
    )
  }, error = function(e) {
    warning("Failed to compute shared controls per treated: ", e$message, call. = FALSE)
    list(mean_shared = NA, median_shared = NA, prop_shared = NA, max_shared = NA)
  })

  return(results)
}

#' Calculate Overlap Statistics from Match Object
#'
#' Extracts a matched table from a match object and then calculates all three
#' overlap metrics.
#'
#' @param mtch A match object from matching procedures.
#'
#' @return A list containing all three sets of overlap statistics.
#'
#' @export
calculate_overlap_statistics_from_match_object <- function(mtch) {
  # Get the full matched table from the match object
  full_matched_table <- tryCatch({
    result_table(mtch, nonzero_weight_only = TRUE)
  }, error = function(e) {
    warning("Failed to get full unit table for overlap statistics: ", e$message, call. = FALSE)
    return(NULL)
  })

  # Then call the function to compute the statistics from the table
  return(calculate_overlap_stats_from_table(full_matched_table))
}


# ============================================================================
# BACKWARD COMPATIBILITY WRAPPER
# ============================================================================

#' Calculate Overlap Statistics (Backward Compatibility)
#'
#' Wrapper function that maintains backward compatibility with the original
#' function name but now returns only the pairwise overlap statistics.
#'
#' @param mtch A match object from matching procedures
#'
#' @return Original pairwise overlap statistics for backward compatibility
#'
#' @export
calculate_overlap_statistics <- function(mtch) {
  all_stats <- calculate_overlap_statistics_from_match_object(mtch)
  return(all_stats$pairwise_overlap)
}
