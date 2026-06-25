"""
Modeling utilities for mangrove canopy height and AGB regression.

Provides a Random Forest training wrapper with RandomizedSearchCV
hyperparameter tuning, regression evaluation metrics, and the standard
diagnostic plots used across this project: 1:1 observed vs predicted
scatter, feature importance, and Pearson correlation matrix.

All randomized procedures use a fixed random_state for reproducibility.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import randint
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

RANDOM_STATE = 42

# Gold-standard RF hyperparameter search space.
RF_PARAM_DIST = {
    'n_estimators': randint(100, 600),
    'max_depth': randint(3, 20),  # None is also a valid manual addition
    'min_samples_split': randint(2, 20),
    'min_samples_leaf': randint(1, 10),
    'max_features': ['sqrt', 'log2', 0.5, 0.8],
    'bootstrap': [True, False],
}


# ============================================================================
# 1. TRAINING
# ============================================================================
def train_rf_random_search(
    X_train,
    y_train,
    param_dist=None,
    n_iter=100,
    cv=5,
    scoring='r2',
    n_jobs=-1,
    random_state=RANDOM_STATE,
    verbose=1,
):
    """
    Train a Random Forest regressor with RandomizedSearchCV tuning.

    Args:
        X_train (array-like): Training features.
        y_train (array-like): Training target.
        param_dist (dict, optional): Hyperparameter distribution.
            Defaults to RF_PARAM_DIST.
        n_iter (int): Number of parameter settings sampled.
        cv (int): Number of cross-validation folds.
        scoring (str): Scoring metric for search.
        n_jobs (int): Parallel jobs (-1 = all cores).
        random_state (int): Random seed for reproducibility.
        verbose (int): Verbosity level passed to RandomizedSearchCV.

    Returns:
        RandomizedSearchCV: Fitted search object. Use `.best_estimator_`
        for the tuned model and `.best_params_` for the winning config.
    """
    param_dist = param_dist or RF_PARAM_DIST

    rf = RandomForestRegressor(random_state=random_state)

    search = RandomizedSearchCV(
        estimator=rf,
        param_distributions=param_dist,
        n_iter=n_iter,
        cv=cv,
        scoring=scoring,
        n_jobs=n_jobs,
        random_state=random_state,
        verbose=verbose,
    )

    search.fit(X_train, y_train)

    print('-' * 30)
    print(f'Best CV {scoring}: {search.best_score_:.4f}')
    print('Best params:')
    for key, val in search.best_params_.items():
        print(f'  {key:<20}: {val}')

    return search


# ============================================================================
# 2. EVALUATION
# ============================================================================
def evaluate_regression(y_true, y_pred, unit='', label='Model'):
    """
    Compute and print standard regression metrics: R2, RMSE, MAE, Bias.

    Args:
        y_true (array-like): Observed values.
        y_pred (array-like): Predicted values.
        unit (str): Unit string for display (e.g. 'm', 'Mg/ha').
        label (str): Label for the printed header.

    Returns:
        dict: {'r2', 'rmse', 'mae', 'bias'}
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    bias = np.mean(y_pred - y_true)

    print(f'--- {label} Regression Metrics (test set) ---')
    print(f'R2  : {r2:.4f}')
    print(f'RMSE: {rmse:.3f} {unit}')
    print(f'MAE : {mae:.3f} {unit}')
    print(f'Bias: {bias:.3f} {unit}')

    return {'r2': r2, 'rmse': rmse, 'mae': mae, 'bias': bias}


def get_feature_importance(model, feature_names):
    """
    Extract feature importance from a fitted tree-based model as a sorted
    DataFrame (descending importance).

    Args:
        model: Fitted estimator with `.feature_importances_`.
        feature_names (list[str]): Names matching the model's input columns.

    Returns:
        pd.DataFrame: Columns ['feature', 'importance'], sorted descending.
    """
    importance_df = pd.DataFrame(
        {'feature': feature_names, 'importance': model.feature_importances_}
    ).sort_values('importance', ascending=False).reset_index(drop=True)

    return importance_df


# ============================================================================
# 3. PLOTTING
# ============================================================================
def plot_scatter_1to1(
    y_true,
    y_pred,
    unit='',
    title='Observed vs Predicted',
    metrics=None,
    ax=None,
    color='#1a73e8',
):
    """
    Plot a 1:1 observed vs predicted scatter plot with regression trendline.

    Args:
        y_true (array-like): Observed values.
        y_pred (array-like): Predicted values.
        unit (str): Unit string for axis labels.
        title (str): Plot title.
        metrics (dict, optional): Output of evaluate_regression(), annotated
            on the plot if provided.
        ax (matplotlib.axes.Axes, optional): Axis to plot on. Creates a new
            figure if None.
        color (str): Point/line color.

    Returns:
        matplotlib.axes.Axes
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    sns.regplot(
        x=y_true, y=y_pred, ax=ax,
        scatter_kws={'alpha': 0.4, 's': 18, 'color': color},
        line_kws={'color': '#e53935', 'linewidth': 1.5},
    )

    lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
    ax.plot(lims, lims, linestyle='--', color='gray', linewidth=1, label='1:1 line')

    ax.set_xlabel(f'Observed ({unit})' if unit else 'Observed')
    ax.set_ylabel(f'Predicted ({unit})' if unit else 'Predicted')
    ax.set_title(title)
    ax.legend(loc='upper left', fontsize=9)

    if metrics:
        text = (
            f"R2 = {metrics['r2']:.3f}\n"
            f"RMSE = {metrics['rmse']:.3f}\n"
            f"MAE = {metrics['mae']:.3f}\n"
            f"Bias = {metrics['bias']:.3f}"
        )
        ax.text(
            0.95, 0.05, text, transform=ax.transAxes,
            fontsize=9, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
        )

    return ax


def plot_feature_importance(importance_df, title='Feature Importance', ax=None, top_n=None):
    """
    Plot a horizontal bar chart of feature importance, most important on top.

    Args:
        importance_df (pd.DataFrame): Output of get_feature_importance().
        title (str): Plot title.
        ax (matplotlib.axes.Axes, optional): Axis to plot on.
        top_n (int, optional): Limit to top N features.

    Returns:
        matplotlib.axes.Axes
    """
    df = importance_df if top_n is None else importance_df.head(top_n)
    df = df.sort_values('importance', ascending=True)

    if ax is None:
        _, ax = plt.subplots(figsize=(6, max(4, 0.35 * len(df))))

    ax.barh(df['feature'], df['importance'], color='#2e7d32')
    ax.set_xlabel('Importance')
    ax.set_title(title)

    return ax


def plot_correlation_matrix(df, columns=None, title='Pearson Correlation Matrix', ax=None):
    """
    Plot an annotated Pearson correlation heatmap for the given columns.

    Args:
        df (pd.DataFrame): Data containing the feature (and optionally
            target) columns.
        columns (list[str], optional): Columns to include. Defaults to all
            numeric columns.
        title (str): Plot title.
        ax (matplotlib.axes.Axes, optional): Axis to plot on.

    Returns:
        matplotlib.axes.Axes
    """
    data = df[columns] if columns is not None else df.select_dtypes(include=np.number)
    corr = data.corr(method='pearson')

    if ax is None:
        _, ax = plt.subplots(figsize=(0.6 * len(corr.columns) + 2, 0.6 * len(corr.columns) + 1))

    sns.heatmap(
        corr, ax=ax, cmap='RdBu_r', vmin=-1, vmax=1, center=0,
        annot=True, fmt='.2f', annot_kws={'size': 8},
        square=True, cbar_kws={'shrink': 0.8},
    )
    ax.set_title(title)

    return ax


def plot_model_diagnostics(
    y_true,
    y_pred,
    importance_df,
    unit='',
    label='Model',
    metrics=None,
    figsize=(12, 5),
):
    """
    Convenience wrapper: side-by-side 1:1 scatter and feature importance
    plot for a single trained model.

    Args:
        y_true (array-like): Observed test values.
        y_pred (array-like): Predicted test values.
        importance_df (pd.DataFrame): Output of get_feature_importance().
        unit (str): Unit string for axis labels.
        label (str): Model label used in titles.
        metrics (dict, optional): Output of evaluate_regression().
        figsize (tuple): Figure size.

    Returns:
        matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    plot_scatter_1to1(
        y_true, y_pred, unit=unit, title=f'{label}: Observed vs Predicted',
        metrics=metrics, ax=axes[0],
    )
    plot_feature_importance(
        importance_df, title=f'{label}: Feature Importance', ax=axes[1],
    )

    fig.tight_layout()
    return fig
